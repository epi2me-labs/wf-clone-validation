#!/usr/bin/env python

"""Report component for displaying information from wf-clone-validation."""
import argparse

from aplanat import base, hist, report
import aplanat.graphics
from aplanat.util import Colors
from bokeh.layouts import layout
import numpy as np
import pandas as pd
from plannotate import annotate
from plannotate import BLAST_hit_details
from plannotate import get_bokeh
import pysam


def run_plannotate(fasta, blast_db, linear=False):
    """Run annotate and create Bokeh plot."""
    with pysam.FastxFile(fasta) as fh:
        seq = next(fh).sequence

    df = annotate(seq, blast_db, linear)
    df = BLAST_hit_details.details(df)
    plot = get_bokeh(df, linear)
    plot.xgrid.grid_line_color = None
    plot.ygrid.grid_line_color = None
    df = clean_results(df)
    return plot, df


def clean_results(df):
    """Clean-up annotation dataframe for display."""
    rename = {
        'Feature': 'Feature', 'uniprot': 'Uniprot ID',
        'db': 'Database',
        'pident': 'Identity', 'abs percmatch': 'Match Length',
        'Description': 'Description',
        'qstart': 'Start Location', 'qend': 'End Location', 'length': 'Length',
        'sframe': 'Strand', 'db': 'Database',
        'qlen': 'qlen'}
    numeric_columns = ['Identity', 'Match Length']
    display_columns = [
        'Feature', 'Uniprot ID',
        'Database',
        'Identity', 'Match Length',
        'Description',
        'Start Location', 'End Location', 'Length',
        'Strand', 'qlen']
    df = df.rename(columns=rename)[display_columns]
    df['Plasmid length'] = df.iloc[0]['qlen']
    df = df.drop(columns='qlen')
    df[numeric_columns] = np.round(df[numeric_columns], 1).astype(str) + "%"
    df.loc[df['Database'] == "infernal", 'Identity'] = "-"
    df.loc[df['Database'] == "infernal", 'Match Length'] = "-"
    df = df.set_index("Feature", drop=True).reset_index()
    return df


def dotplot_assembly(assembly_maf):
    """Dotplot the assembly."""
    records = list()
    with open(assembly_maf) as maf:
        while True:
            line = maf.readline()
            print(line)
            if line.startswith('#'):
                continue
            elif line.startswith('a'):
                r1 = maf.readline().split()[1:5]
                r2 = maf.readline().split()[1:5]
                maf.readline()
                records.append(r1 + r2)
            elif line == "":
                break
            else:
                print(line)
                raise IOError("Cannot read alignment file")

    names = ['ref', 'rstart', 'rlen', 'rorient',
             'query', 'qstart', 'qlen', 'qorient']
    df = pd.DataFrame(records, columns=names)
    df = df.loc[df['qorient'] == '+']
    for f in ['qstart', 'qlen', 'rstart', 'rlen']:
        df[f] = df[f].astype(int)
    df['qend'] = df['qstart'] + df['qlen']
    df['rend'] = df['rstart'] + df['rlen']

    df['qend'] = df['qend'].astype(int)
    df['rend'] = df['rend'].astype(int)
    dotplot = base.simple(
        [], [],
        xlim=(0, max(df['rend'])),
        ylim=(0, max(df['qend'])),
        width=440,
        height=400,
        x_axis_label='position',
        y_axis_label='position',
        title=f"Dot plot for: {assembly_maf}"
    )
    dotplot.segment(df['rstart'], df['qstart'], df['rend'], df['qend'])
    return dotplot


def per_assembly(assemblies_file, database, mafs):
    """Build_per_assembly_tuples."""
    tup = []

    for maf in mafs:
        dotplot = dotplot_assembly(maf)
        filename = str(maf)[:-4]
        plot, annotations = run_plannotate(filename, database)
        df = pd.read_csv(assemblies_file, sep='\t', index_col=0)
        num_seqs = df.loc[str(filename), 'num_seqs']
        min_len = df.loc[str(filename), 'min_len']
        avg_len = df.loc[str(filename), 'avg_len']
        max_len = df.loc[str(filename), 'max_len']
        tup.append((filename[:-12:], num_seqs, min_len,
                    avg_len,
                    max_len, plot,
                    dotplot, annotations))

    return(tup)


def plot_read_length_distribution(
            fname, nseqs, nbases, minl, maxl, data):
    """Plot_read_length_distribution."""
    kilobases = int(nbases//1000)
    plot = hist.histogram(
        [data['read_length'].tolist()],
        # bins=1000,
        height=300,
        width=400,
        xlim=(0, max(data['read_length']) + 200),
        colors=[Colors.light_cornflower_blue],
        x_axis_label='Length',
        y_axis_label='Count',
        title=f"{nseqs} reads, {kilobases} Kb total, {minl} min, {maxl} max"
    )
    return plot


def plot_qscore_distribution(fname, mean, data):
    """Plot_qscore_distribution."""
    plot = hist.histogram(
        [data['mean_quality'].tolist()],
        # bins=600,
        height=300,
        width=400,
        xlim=(0, 30),
        colors=[Colors.light_cornflower_blue],
        x_axis_label='Mean Quality',
        y_axis_label='Count',
        title=f"Mean Q-score: {mean}"
    )

    return plot


def build_samples_panel(summary_file, reads_file):
    """Build_length_tab."""
    summary_df = pd.read_csv(summary_file, sep="\t")
    reads_df = pd.read_csv(reads_file, sep="\t")
    qc_plots_dic = {}
    filenames = set(reads_df['filename'].tolist())
    for fname in filenames:
        fname_read = reads_df.loc[reads_df['filename'] == fname]
        fname_summary = summary_df.loc[summary_df['filename'] == fname]
        length_plot = plot_read_length_distribution(
            fname, int(fname_summary['n_seqs']),
            int(fname_summary['n_bases']),
            int(fname_summary['min_length']),
            int(fname_summary['max_length']), fname_read)
        qual_plot = plot_qscore_distribution(
            fname, float(fname_summary['mean_quality']), fname_read)
        file_name = str(fname)[:-18:]
        qc_plots_dic[file_name] = (length_plot, qual_plot)

    return qc_plots_dic


def tidyup_status_file(status_sheet, mafs, sample_sheet):
    """Tidy up the sample status file."""
    status_sheet = open(status_sheet).read().splitlines()
    pass_fail = {}
    # mark samples as passed or failed
    for i in status_sheet[1::]:
        i = i.split('/')
        pass_fail[i[-1]] = 'fail'
    barcode_dic = {}
    # if there is a sample sheet translated barcode>sample otherwise pass
    try:
        sample_sheet = open(sample_sheet).read().splitlines()
        for i in sample_sheet[1::]:
            i = i.split(',')
            barcode_dic[i[0]] = i[1]
        translated = {}
        for k, v in pass_fail.items():
            sample_name = barcode_dic[k]
            pass_or_fail = pass_fail[k]
            translated[sample_name] = pass_or_fail
    except Exception:
        translated = pass_fail

    # change status to pass for any samples that made it to end of pipeline
    for maf in mafs:
        name = str(maf)[:-16]
        translated[name] = 'pass'
    status_df = pd.DataFrame(translated.items(),
                             columns=['Sample', 'pass/fail'])
    status_df.to_csv('sample_status.txt', index=False)
    return(status_df)


def output_feature_table(data):
    """Build feature table text file."""
    df = data[0][7]
    sample_column = data[0][0]
    df.insert(0, 'Sample_name', sample_column)
    df.to_csv('feature_table.txt', mode='a', header=True, index=False)
    for i in data[1:]:
        df = i[7]
        sample_column = i[0]
        df.insert(0, 'Sample_name', sample_column)
        df.to_csv('feature_table.txt', mode='a', header=False, index=False)


def main():
    """Entry point to create a wf-clone-validation report."""
    parser = argparse.ArgumentParser(
        'Clone Validation QC report',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False)
    parser.add_argument(
        "--assembly_summary",
        required=True
    )
    parser.add_argument(
        "--assembly_mafs",
        nargs='*',
        required=True
    )
    parser.add_argument(
        "--reads_summary",
        required=True
    )
    parser.add_argument(
        "--fastq_summary",
        required=True
    )
    parser.add_argument(
        "--consensus",
        nargs='+',
        required=True
    )
    parser.add_argument(
        "--revision", default='unknown',
        help="revision number")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit number")
    parser.add_argument(
        "--database", default='unknown',
        help="database to use, directory containing BLAST et. al. files.")
    parser.add_argument(
        "--status_sheet", default='unknown',
        help="status file")
    parser.add_argument(
        "--sample_sheet",
        help="status file")
    args = parser.parse_args()
    report_doc = report.WFReport(
        "Clone Validation Report", "wf-clone-validation",
        revision=args.revision, commit=args.commit)
    section = report_doc.add_section()
    section.markdown("Results generated through the wf-clone-validation "
                     "nextflow workflow provided by Oxford Nanopore "
                     "Technologies")
    section = report_doc.add_section()
    section.markdown("""
## Assemblies
For each assembly read length statistics are displayed a [pLannotate
plot](https://github.com/barricklab/pLannotate), a plot of the quality, a
dotplot generated by self-alignment.  These dotplots should be a near-perfect
diagonal line with little of no off-diagonal elements.  The feature table
provides descriptions of the annotated sequence.

Unfilled features on the plannotate plots are incomplete features; the sequence
match in the plasmid covers less than 95% of the full length of the feature in
the database. These elements may be leftover fragments from earlier cloning
steps used to create a plasmid. If they include only a small fraction of the
feature, they likely do not still have the annotated function. However, even
small feature fragments may affect plasmid function if they result in cryptic
gene expression or are inadvertently combined with other elements during later
cloning steps.

The Plasmid annotation plot and feature table are produced using
[Plannotate](http://plannotate.barricklab.org/).
""")

    assembly = args.assembly_summary
    mafs = args.assembly_mafs
    fastq_summ = args.fastq_summary
    reads_summ = args.reads_summary
    database = args.database
    alldata = per_assembly(assembly, database, mafs)
    output_feature_table(alldata)
    summary_stats_dic = build_samples_panel(fastq_summ, reads_summ)

    for i in alldata:
        section = report_doc.add_section()
        section.markdown('### Sample: {}'.format(str(i[0])))
        exec_summary = aplanat.graphics.InfoGraphItems()
        exec_summary.append("No. Seqs",
                            str(i[1]),
                            "bars", '')
        exec_summary.append("Min Length",
                            str(i[2]),
                            "align-left", '')
        exec_summary.append("Average Length",
                            str(int(i[3])),
                            "align-center", '')
        exec_summary.append("Max Length",
                            str(i[4]),
                            'align-right')
        exec_plot = aplanat.graphics.infographic(
                    exec_summary.values(), ncols=4)
        section.plot(exec_plot, key="exec-plot"+str(i[1]))
        dotplot = [i[6]]
        lengthplot = summary_stats_dic[str(i[0])][0]
        qstatplot = summary_stats_dic[str(i[0])][1]
        plasmidplot = [i[5]]
        section.plot(
            layout(
                [[lengthplot, qstatplot], [dotplot, plasmidplot]],
                sizing_mode='scale_width'))
        section.table((i[7].drop(columns=['Plasmid length'])),
                      index=False, key="table"+str(i[1]))
    section = report_doc.add_section()
    status_sheet = args.status_sheet
    sample_sheet = args.sample_sheet
    pass_fail = tidyup_status_file(status_sheet, mafs, sample_sheet)
    section.markdown("##Sample status")
    section.table(pass_fail, index=False)

    report_doc.write('wf-clone-validation-report.html')


if __name__ == "__main__":
    main()
