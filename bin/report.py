#!/usr/bin/env python

"""Report component for displaying information from wf-clone-validation."""
import argparse
import json
import os

from aplanat import bars, base, report
from aplanat import json_item
from aplanat.components import fastcat
import aplanat.graphics
from aplanat.util import Colors
from bokeh.layouts import layout
from bokeh.models import Panel, Tabs
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


def per_assembly(assemblies_file, database, sample_files, item):
    """Build_per_assembly_tuples."""
    dotplot = dotplot_assembly(sample_files['maf'])
    plot, annotations = run_plannotate(sample_files['fasta'], database)
    df = pd.read_csv(assemblies_file, sep='\t', index_col=0)
    filename = str(item) + '.final.fasta'
    num_seqs = df.loc[filename, 'num_seqs']
    min_len = df.loc[filename, 'min_len']
    avg_len = df.loc[filename, 'avg_len']
    max_len = df.loc[filename, 'max_len']
    tup = {'sample_name': item,
           'num_seqs': num_seqs,
           'min_len': min_len,
           'avg_len': avg_len,
           'max_len': max_len,
           'plot': plot,
           'dotplot': dotplot,
           'annotations': annotations}
    return(tup)


def tidyup_status_file(status_sheet):
    """Tidy up the sample status file."""
    sample_status = pd.read_csv(status_sheet[0], header=None)
    unique_samples = sample_status[0].unique()
    pass_fail_dic = {}
    for sample in unique_samples:
        pass_fail_dic[sample] = 'Pass'
    filter_pass = sample_status[sample_status[1] != 'Pass']
    failures = dict(zip(filter_pass[0], filter_pass[1]))
    all_sample_names = unique_samples.tolist()
    all_sample_names.sort()
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
    return(status_df, passed_list, all_sample_names)


def output_feature_table(data):
    """Build feature table text file or if no data output empty file."""
    if data:
        df = data[0]['annotations']
        sample_column = data[0]['sample_name']
        df.insert(0, 'Sample_name', sample_column)
        df.to_csv('feature_table.txt', mode='a', header=True, index=False)
        for sample in data[1:]:
            df = sample['annotations']
            sample_column = sample['sample_name']
            df.insert(0, 'Sample_name', sample_column)
            df.to_csv('feature_table.txt', mode='a', header=False, index=False)
    else:
        # If no samples passed create empty feature_table file
        feature_file = open("feature_table.txt", "w")
        feature_file.write(
"""Sample_name,Feature,Uniprot ID,Database,Identity,Match Length,Description,Start Location,End Location,Length,Strand,Plasmid length""") # noqa
        feature_file.close()


def pair_samples_with_mafs(sample_names):
    """Match Assembly sequences with mafs."""
    fasta_mafs = {}
    for sample_name in sample_names:
        fasta = 'assemblies/' + sample_name + '.final.fasta'
        maf = 'assembly_maf/' + sample_name + '.final.fasta.maf'
        if os.path.exists(fasta) and os.path.exists(maf):
            fasta_mafs[sample_name] = {'fasta': fasta,
                                       'maf': maf}
        else:
            print("Missing data required for report: " + sample_name)
    return(fasta_mafs)


def read_files(summaries, sep='\t'):
    """Read a set of files and join to single dataframe."""
    dfs = list()
    for fname in sorted(summaries):
        dfs.append(pd.read_csv(fname, sep=sep))
    return pd.concat(dfs)


def fastcat_report_tab(file_name, tab_name):
    """Read fastcat dataframe and create a tab with qual and len plots."""
    df = pd.read_csv(file_name, sep='\t')
    depth = len(df.index)
    min_length = df["read_length"].min()
    max_length = df["read_length"].max()
    lengthplot = fastcat.read_length_plot(
                df,
                min_len=min_length,
                max_len=max_length)
    qstatplot = fastcat.read_quality_plot(df)
    exec_summary = aplanat.graphics.InfoGraphItems()
    exec_summary.append('Total No. samples',
                        str(depth),
                        "bars", '')
    exec_plot = aplanat.graphics.infographic(
        exec_summary.values(), ncols=1)
    tab = Panel(child=layout(
        [[exec_plot], [lengthplot, qstatplot]],
        aspect_ratio="auto",
        sizing_mode='stretch_width'),
        title=tab_name)
    return tab


def create_fastcat_dic(sample_names, raw, hostfilt, downsampled):
    """Create dictionary using sample names and fastcat files available."""
    per_sample_dic = {}
    lists = {'raw': raw, 'hostfilt': hostfilt, 'downsampled': downsampled}
    for sample in sample_names:
        new_dic = {}
        item_search = '/' + sample + '.'
        for list_name, fc_list in lists.items():
            #  find index of item that contains the sample name as a substring
            item_index = [i for i, s in enumerate(fc_list) if item_search in s]
            if item_index:
                indice = item_index[0]
                new_dic[list_name] = fc_list[indice]
            else:
                pass
        per_sample_dic[sample] = new_dic
    return per_sample_dic


def exec_summary_plot(tup_dic):
    """Create and return the infographic summary."""
    exec_summary = aplanat.graphics.InfoGraphItems()
    exec_summary.append(
        "No. Seqs",
        str(tup_dic['num_seqs']),
        "bars", '')
    exec_summary.append(
        "Min Length",
        str(tup_dic['min_len']),
        "align-left", '')
    exec_summary.append(
        "Average Length",
        str(int(tup_dic['avg_len'])),
        "align-center", '')
    exec_summary.append(
        "Max Length",
        str(tup_dic['max_len']),
        'align-right')
    exec_plot = aplanat.graphics.infographic(
        exec_summary.values(), ncols=4)
    return(exec_plot)


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
        "--downsampled_stats",
        nargs='+',
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
    parser.add_argument(
        "--per_barcode_stats", nargs='+',
        help="fastcat stats file for each sample before filtering")
    parser.add_argument(
        "--host_filter_stats", nargs='+',
        help="fastcat stats file after host filtering")
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
    seq_summary = read_files(args.per_barcode_stats)
    section = report_doc.add_section()
    section.markdown("###Number of reads per barcode")
    barcode_counts = (
        pd.DataFrame(seq_summary['sample_name'].value_counts())
        .sort_index()
        .reset_index()
        .rename(
            columns={'index': 'sample', 'sample_name': 'count'})
    )
    bc_counts = bars.simple_bar(
        barcode_counts['sample'].astype(str), barcode_counts['count'],
        colors=[Colors.cerulean]*len(barcode_counts),
        title=(
            'Number of reads per barcode'),
        plot_width=None
    )
    bc_counts.xaxis.major_label_orientation = 3.14/2
    section.plot(
        layout(
            [[bc_counts]],
            sizing_mode="stretch_width"))
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
    status = args.status
    pass_fail = tidyup_status_file(status)
    sample_data = pair_samples_with_mafs(pass_fail[1])
    passed_samples = pass_fail[1]
    host_ref_stats = args.host_filter_stats
    downsampled_stats = args.downsampled_stats
    database = args.database
    initial_stats = args.per_barcode_stats
    if ('host_filter_stats/OPTIONAL_FILE' in host_ref_stats):
        host_filt = host_ref_stats.remove('host_filter_stats/OPTIONAL_FILE')
        if host_filt is None:
            host_filt = []
    else:
        host_filt = host_ref_stats
    if ('host_filter_stats/OPTIONAL_FILE' in downsampled_stats):
        summary_stats = downsampled_stats.remove(
                        'downsampled_stats/OPTIONAL_FILE')
        if summary_stats is None:
            summary_stats = []
    else:
        summary_stats = downsampled_stats
    sample_names = pass_fail[2]
    final_samples = []
    fast_cat_dic = create_fastcat_dic(sample_names, initial_stats,
                                      host_filt, summary_stats)
    json_file = open("plannotate.json", "a")
    plannotate_collection = {}
    for item in sample_names:
        if item in passed_samples:
            sample_files = sample_data[item]
            tup_dic = per_assembly(assembly, database, sample_files, item)
            final_samples.append(tup_dic)
            section = report_doc.add_section()
            section.markdown('### Sample: {}'.format(str(item)))
            infographic_plot = exec_summary_plot(tup_dic)
            section.plot(infographic_plot, key="exec-plot"+(str(item)))
            fast_cat_tabs = fast_cat_dic[item]
            alltabs = []
            for key, value in fast_cat_tabs.items():
                alltabs.append(fastcat_report_tab(value, key))
            cover_panel = Tabs(tabs=alltabs)
            dotplot = [tup_dic['dotplot']]
            plasmidplot = [tup_dic['plot']]
            section.plot(layout(
                            [cover_panel,
                                [dotplot, plasmidplot]],
                            sizing_mode='scale_width'))
            section.table((tup_dic['annotations']).drop(
                columns=['Plasmid length']),
                index=False, key="table"+str(item))
            plasmid_len = tup_dic['annotations']['Plasmid length'][0]
            plannotate_dic = {"barcode": item, "reflen": plasmid_len}
            feature_dic = tup_dic['annotations'].drop(
                          ['Plasmid length'], axis=1)
            features = feature_dic.to_dict('records')
            output_json = json_item(tup_dic['plot'])
            plannotate_dic = {"reflen": float(plasmid_len),
                              "features": features,
                              "plot": output_json['doc']}
            plannotate_collection[item] = plannotate_dic
        else:
            section.markdown('### Sample Failed: {}'.format(str(item)))
            fast_cat_tabs = fast_cat_dic[item]
            alltabs = []
            for key, value in fast_cat_tabs.items():
                alltabs.append(fastcat_report_tab(value, key))
            cover_panel = Tabs(tabs=alltabs)
            section.plot(cover_panel)
    json_object = json.dumps(plannotate_collection, indent=4)
    json_file.write(json_object)
    json_file.close()
    output_feature_table(final_samples)
    section = report_doc.add_section()
    status = args.status
    pass_fail = tidyup_status_file(status)
    section.markdown("##Sample status")
    section.table(pass_fail[0], index=False)

    report_doc.write('wf-clone-validation-report.html')


if __name__ == "__main__":
    main()
