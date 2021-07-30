#!/usr/bin/env python

"""Report component for displaying information from wf-clone-validation."""
import argparse
import os

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


def per_assembly(assemblies_file, database, sample_data):
    """Build_per_assembly_tuples."""
    tup = []
    for sample in sample_data:
        dotplot = dotplot_assembly(sample['maf'])
        plot, annotations = run_plannotate(sample['fasta'], database)
        df = pd.read_csv(assemblies_file, sep='\t', index_col=0)
        filename = str(sample['sample_name']) + '.final.fasta'
        num_seqs = df.loc[filename, 'num_seqs']
        min_len = df.loc[filename, 'min_len']
        avg_len = df.loc[filename, 'avg_len']
        max_len = df.loc[filename, 'max_len']
        tup.append({'sample_name': sample['sample_name'],
                    'num_seqs': num_seqs,
                    'min_len': min_len,
                    'avg_len': avg_len,
                    'max_len': max_len,
                    'plot': plot,
                    'dotplot': dotplot,
                    'annotations': annotations})
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


def tidyup_status_file(status_sheet):
    """Tidy up the sample status file."""
    sample_status = pd.read_csv(status_sheet[0], header=None)
    unique_samples = sample_status[0].unique()
    pass_fail_dic = {}
    for sample in unique_samples:
        pass_fail_dic[sample] = 'Pass'
    filter_pass = sample_status[sample_status[1] != 'Pass']
    failures = dict(zip(filter_pass[0], filter_pass[1]))
    passed_list = unique_samples.tolist()
    for k, v in failures.items():
        pass_fail_dic[k] = v
        passed_list.remove(k)
    passed_list.sort()
    status_df = pd.DataFrame(pass_fail_dic.items(),
                             columns=['Sample', 'pass/failed reason'])
    sort_df = status_df['Sample'].astype(str).argsort()
    status_df = status_df.iloc[sort_df]
    status_df.to_csv('sample_status.txt', index=False)
    return(status_df, passed_list)


def output_feature_table(data):
    """Build feature table text file."""
    df = data[0]['annotations']
    sample_column = data[0]['sample_name']
    df.insert(0, 'Sample_name', sample_column)
    df.to_csv('feature_table.txt', mode='a', header=True, index=False)
    for sample in data[1:]:
        df = sample['annotations']
        sample_column = sample['sample_name']
        df.insert(0, 'Sample_name', sample_column)
        df.to_csv('feature_table.txt', mode='a', header=False, index=False)


def pair_samples_with_mafs(sample_names):
    """Match Assembly sequences with mafs."""
    fasta_mafs = []
    for sample_name in sample_names:
        fasta = 'assemblies/' + sample_name + '.final.fasta'
        maf = 'assembly_maf/' + sample_name + '.final.fasta.maf'
        if os.path.exists(fasta) and os.path.exists(maf):
            fasta_mafs.append({'sample_name': sample_name,
                               'fasta': fasta,
                               'maf': maf})
        else:
            print("Missing data required for report: " + sample_name)
    return(fasta_mafs)


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
        "--status", nargs='+',
        help="status")
    args = parser.parse_args()
    report_doc = report.WFReport(
        "Clone Validation Report",
        "wf-clone-validation",
        revision=args.revision,
        commit=args.commit)
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
    status = args.status
    pass_fail = tidyup_status_file(status)
    sample_data = pair_samples_with_mafs(pass_fail[1])
    fastq_summ = args.fastq_summary
    reads_summ = args.reads_summary
    database = args.database
    if (mafs[0] != 'assembly_maf/OPTIONAL_FILE'):
        alldata = per_assembly(assembly, database, sample_data)
        output_feature_table(alldata)
        summary_stats_dic = build_samples_panel(fastq_summ, reads_summ)
        for i in alldata:
            section = report_doc.add_section()
            section.markdown('### Sample: {}'.format(i['sample_name']))
            exec_summary = aplanat.graphics.InfoGraphItems()
            exec_summary.append(
                "No. Seqs",
                str(i['num_seqs']),
                "bars", '')
            exec_summary.append(
                "Min Length",
                str(i['min_len']),
                "align-left", '')
            exec_summary.append(
                "Average Length",
                str(int(i['avg_len'])),
                "align-center", '')
            exec_summary.append(
                "Max Length",
                str(i['max_len']),
                'align-right')
            exec_plot = aplanat.graphics.infographic(
                exec_summary.values(), ncols=4)
            section.plot(exec_plot, key="exec-plot"+(i['sample_name']))
            dotplot = [i['dotplot']]
            lengthplot = summary_stats_dic[str(i['sample_name'])][0]
            qstatplot = summary_stats_dic[i['sample_name']][1]
            plasmidplot = [i['plot']]
            section.plot(
                layout(
                    [[lengthplot, qstatplot], [dotplot, plasmidplot]],
                    sizing_mode='scale_width'))
            section.table((i['annotations'].drop(columns=['Plasmid length'])),
                          index=False, key="table"+(str(i['sample_name'])))
    else:
        open("feature_table.txt", "w")

    section = report_doc.add_section()
    status = args.status
    pass_fail = tidyup_status_file(status)
    section.markdown("##Sample status")
    section.table(pass_fail[0], index=False)

    report_doc.write('wf-clone-validation-report.html')


if __name__ == "__main__":
    main()
