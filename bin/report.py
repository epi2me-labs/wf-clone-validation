#!/usr/bin/env python

"""Report component for displaying information from wf-clone-validation."""
import argparse
import json
import os


from aplanat import bars, report
from aplanat.components import fastcat
from aplanat.components import simple as scomponents
import aplanat.graphics
from aplanat.util import Colors
from bokeh.layouts import layout
from bokeh.models import Panel, Tabs
import pandas as pd
from plannotate import get_bokeh


def tidyup_status_file(status_sheet):
    """Tidy up the sample status file."""
    sample_status = pd.read_csv(status_sheet[0], header=None)
    unique_samples = sample_status[0].unique()
    pass_fail_dic = {}
    for sample in unique_samples:
        pass_fail_dic[sample] = 'Completed successfully'
    filter_pass = sample_status[sample_status[1] != 'Completed successfully']
    failures = dict(zip(filter_pass[0], filter_pass[1]))
    all_sample_names = unique_samples.tolist()
    all_sample_names.sort()
    passed_list = unique_samples.tolist()
    for k, v in failures.items():
        pass_fail_dic[k] = v
        if v != 'Completed but failed to reconcile':
            passed_list.remove(k)
    passed_list.sort()
    status_df = pd.DataFrame(
        pass_fail_dic.items(), columns=['Sample', 'pass/failed reason'])
    sort_df = status_df['Sample'].astype(str).argsort()
    status_df = status_df.iloc[sort_df]
    return(status_df, passed_list, all_sample_names, pass_fail_dic)


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
    exec_summary.append(
        'No. reads',
        str(depth),
        "bars", '')
    exec_plot = aplanat.graphics.infographic(
        exec_summary.values(), ncols=1)
    tab = Panel(child=layout(
        [[exec_plot], [lengthplot], [qstatplot]],
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


def main():
    """Entry point to create a wf-clone-validation report."""
    parser = argparse.ArgumentParser(
        'Clone Validation QC report',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False)
    parser.add_argument(
        "--downsampled_stats",
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
        "--status", nargs='+',
        help="status")
    parser.add_argument(
        "--per_barcode_stats", nargs='+',
        help="fastcat stats file for each sample before filtering")
    parser.add_argument(
        "--host_filter_stats", nargs='+',
        help="fastcat stats file after host filtering")
    parser.add_argument(
        "--params", default=None,
        help="A csv containing the parameter key/values")
    parser.add_argument(
        "--versions",
        help="directory contained CSVs containing name,version.")
    parser.add_argument(
        "--plannotate_json",
        help="Plannotate Json.")
    parser.add_argument(
        "--inserts_json",
        help="inserts Json.")
    parser.add_argument(
        "--report_name",
        help="report name")
    args = parser.parse_args()
    report_doc = report.WFReport(
        "Clone Validation Report",
        "wf-clone-validation",
        revision=args.revision,
        commit=args.commit)

    # summary section
    section = report_doc.add_section()
    section.markdown("### Summary")
    seq_summary = read_files(args.per_barcode_stats)
    section = report_doc.add_section()
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

    # We defer this until processing through all the samples in the loop below
    summary_placeholder = report_doc.add_section(key='stats_table')
    pass_fail = tidyup_status_file(args.status)
    # find inserts and put in dir
    current_directory = os.getcwd()
    make_directory = os.path.join(current_directory, r'inserts')
    if not os.path.exists(make_directory):
        os.makedirs(make_directory)
    if os.stat(args.inserts_json).st_size != 0:
        inserts = json.load(open(args.inserts_json))
        insert_placeholder = report_doc.add_section(key='insert_table')
        seq_segment = inserts['bed_df']
        seq_segment = pd.read_json(inserts['bed_df'])
        insert_placeholder.markdown("""
### Insert sequences
This table shows which primers were found in the consensus sequence
of each sample and where the inserts were found.
""")
        section = report_doc.add_section()
        section.markdown("""
### Multiple sequence alignment
This section shows the inserts aligned with each other or a reference
sequence if provided.
""")
        section.markdown("<pre>" + os.linesep.join(inserts['msa']) + "</pre>")

    # Per sample details
    section = report_doc.add_section()
    section.markdown("""
### Assemblies
For each assembly read length statistics are displayed a [pLannotate
plot](https://github.com/barricklab/pLannotate), a plot of the quality.
The feature table provides descriptions of the annotated sequence.

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
    passed_samples = pass_fail[1]
    host_ref_stats = args.host_filter_stats
    downsampled_stats = args.downsampled_stats
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
    fast_cat_dic = create_fastcat_dic(
        sample_names, initial_stats, host_filt, summary_stats)
    plannotate = json.load(open(args.plannotate_json))
    sample_stats = []

    # stats graphs and plannotate where appropriate
    for item in sample_names:
        section = report_doc.add_section()
        section.markdown('### Sample: {}'.format(str(item)))
        if item in passed_samples:
            section.markdown('*{}*'.format(pass_fail[3][item]))
            tup_dic = plannotate[item]
            annotations = pd.read_json(tup_dic['annotations'])
            plot = get_bokeh(pd.read_json(tup_dic['plot']), False)
            fast_cat_tabs = fast_cat_dic[item]
            alltabs = []
            for key, value in fast_cat_tabs.items():
                alltabs.append(fastcat_report_tab(value, key))
            cover_panel = Tabs(tabs=alltabs)
            plasmidplot = [plot]
            stats_table = [item] + [int(tup_dic['seq_len'])]
            sample_stats.append(stats_table)
            section.plot(layout(
                            [[cover_panel, plasmidplot]],
                            sizing_mode='scale_width'))
            section.table(annotations.drop(
                columns=['Plasmid length']),
                index=False, key="table"+str(item))
        else:
            section.markdown('*{}*'.format(pass_fail[3][item]))
            fast_cat_tabs = fast_cat_dic[item]
            alltabs = []
            for key, value in fast_cat_tabs.items():
                alltabs.append(fastcat_report_tab(value, key))
            cover_panel = Tabs(tabs=alltabs)
            section.plot(cover_panel)
            stats_table = [item] + ['N/A']
            sample_stats.append(stats_table)

    # high level sample status table
    summary_placeholder.markdown("### Sample status")
    stats_df = pd.DataFrame(
        sample_stats, columns=["Sample", "Length"])
    merged_inner = pd.merge(pass_fail[0], stats_df)
    summary_placeholder.table(merged_inner, index=False, key='stats_table')
    merged_inner.to_csv('sample_status.txt', index=False)

    def insert_len(start, end, length):
        if end >= start:
            insert_length = end - start
        else:
            insert_length = (length - start) + end
        return insert_length

    if os.stat(args.inserts_json).st_size != 0:
        test_df = seq_segment.merge(stats_df, how='left')
        test_df['Insert length'] = list(
            map(
                insert_len,
                test_df['start'],
                test_df['end'],
                test_df['Length']))
        test_df = test_df.drop(columns='Length')
        insert_placeholder.table(test_df, index=False, key='insert_table')
    # Versions and params
    report_doc.add_section(
        section=scomponents.version_table(args.versions))
    report_doc.add_section(
        section=scomponents.params_table(args.params))
    report_doc.write(args.report_name)


if __name__ == "__main__":
    main()
