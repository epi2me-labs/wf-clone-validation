#!/usr/bin/env python
"""Create workflow report."""
from itertools import count, takewhile
import json
import os

from dominate.tags import p, pre
from dominate.util import raw
from ezcharts.components import bokehchart
from ezcharts.components.fastcat import draw_all_plots
from ezcharts.components.reports import labs
from ezcharts.layout.snippets import Stats, Tabs
from ezcharts.layout.snippets.table import DataTable
from ezcharts.plots.__init__ import BokehPlot
from ezcharts.plots.util import read_files
from ezcharts.util import get_named_logger  # noqa: ABS101
import pandas as pd
from workflow_glue.bokeh_plot import get_bokeh
from workflow_glue.report_utils import report_utils

from .util import wf_parser  # noqa: ABS101

THEME = 'epi2melabs'


def sliced(seq, n):
    """Slice up a string."""
    iterator = takewhile(len, (seq[i: i + n] for i in count(0, n)))
    slices = list(iterator)
    return slices


def format_msa(inserts):
    """Format the MSA."""
    collect_sliced = []
    name_list = []
    reorg_slices = []
    for i in list(inserts['msa']):
        name_list.append(i[0:10])
        collect_sliced.append(sliced(i[10:], 80))
    for slice_index in range(len(collect_sliced[0])):
        for index, value in enumerate(collect_sliced):
            reorg_slices.append(name_list[index] + value[slice_index])
        reorg_slices.append('')
    return reorg_slices


def button_format(reason):
    """Format button, blue for success, red for failure."""
    if reason == "Completed successfully":
        return """<span class="badge bg-primary">""" + reason + """</span>"""
    else:
        return """<span class="badge bg-danger">""" + reason + """</span>"""


def main(args):
    """Run the entry point."""
    logger = get_named_logger("Report")
    report = labs.LabsReport(
        "Clone validation report", "wf-clone-validation-new",
        args.params, args.versions)
    plannotate_annotations = json.load(open(args.plannotate_json)).keys()
    # Get sample info from status file
    passed_samples, sample_names, sample_status_dic = \
        report_utils.tidyup_status_file(
            args.status, plannotate_annotations)
    # Sample status table
    with report.add_section("Sample status", "Sample status"):
        p("""
This table gives the status of each sample, either completed \
successfully or the reason the workflow failed. \
If applicable the mean quality for the whole \
plasmid assembly has been provided, derived by Medaka \
from aligning the reads to the consensus.
""")
        lengths = json.load(open(args.lengths))
        sample_stats = []
        for item in sample_names:
            try:
                stats_table = [item] + [int(lengths[item]['reflen'])]
            except KeyError:
                stats_table = [item] + ['N/A']
            sample_stats.append(stats_table)
        stats_df = pd.DataFrame(
            sample_stats, columns=["Sample", "Length"])
        status_df = pd.DataFrame(
            sample_status_dic.items(), columns=['Sample', 'pass/failed reason'])
        sort_df = status_df['Sample'].astype(str).argsort()
        status_df = status_df.iloc[sort_df]
        merged = pd.merge(status_df, stats_df)
        merged.to_csv('sample_status.txt', index=False)
        if ('assembly_quality/OPTIONAL_FILE' not in args.assembly_quality):
            qc_df = read_files(args.assembly_quality)[
                ['sample_name', 'mean_quality']]
            qc_df = qc_df.rename(columns={
                'sample_name': 'Sample', 'mean_quality': 'Mean Quality'})
            qc_df = qc_df.reset_index(drop=True)
            merged = pd.merge(merged, qc_df, how="outer")
            merged.fillna('N/A', inplace=True)
            merged['pass/failed reason'] = merged.apply(
                lambda merged: button_format(merged['pass/failed reason']),
                axis=1)
        DataTable.from_pandas(merged, use_index=False)
    with report.add_section("Plannotate", "Plannotate"):
        raw("""The Plasmid annotation plot and feature table are produced using \
<a href="http://plannotate.barricklab.org/">pLannotate</a>""")
        p("""
A pLannotate plot is shown for each assembly.
A feature table provides descriptions of the annotated sequence.

Unfilled features on the plannotate plots are incomplete features; the sequence
match in the plasmid covers less than 95% of the full length of the feature in
the database. These elements may be leftover fragments from earlier cloning
steps used to create a plasmid. If they include only a small fraction of the
feature, they likely do not still have the annotated function. However, even
small feature fragments may affect plasmid function if they result in cryptic
gene expression or are inadvertently combined with other elements during later
cloning steps.
""")
        # Load and plot per sample plannotate json
        plannotate = json.load(open(args.plannotate_json))
        tabs = Tabs()
        with tabs.add_dropdown_menu():
            for item in passed_samples:
                plan_item = plannotate[item]
                annotations = pd.read_json(plan_item['annotations'])
                annotations = annotations.drop("Plasmid length", axis='columns')
                annotations.loc[annotations['Strand'] == -1, 'Strand'] = "-"
                annotations.loc[annotations['Strand'] == 1, 'Strand'] = "+"
                with tabs.add_dropdown_tab(str(item)):
                    bk_plot = BokehPlot()
                    bk_plot._fig = get_bokeh(pd.read_json(plan_item['plot']))
                    bk_plot._fig.xgrid.grid_line_color = None
                    bk_plot._fig.ygrid.grid_line_color = None
                    bokehchart.BokehChart(
                        bk_plot, 'epi2melabs',
                        width='100%', height='100%')
                    DataTable.from_pandas(annotations, use_index=False)
    # Per barcode read count plot
    report_utils.read_count_barplot(args.per_barcode_stats, report)
    with report.add_section("Read stats", "Read stats"):
        p("""
For each assembly, read length statistics and plots of quality \
before and after downsampling, and after host filtering \
if a host reference was provided.
""")
        fastcat_dic = report_utils.create_fastcat_dic(
            sample_names, args.per_barcode_stats,
            args.host_filter_stats,  args.downsampled_stats)
        tabs = Tabs()
        with tabs.add_dropdown_menu():
            for item in sample_names:
                with tabs.add_dropdown_tab(str(item)):
                    p(sample_status_dic[item])
                    internal_tabs = Tabs()
                    # For each sample, plot fastcat stats
                    fastcat_tabs = fastcat_dic[item]
                    for key, value in fastcat_tabs.items():
                        with internal_tabs.add_tab(key):
                            seq_summary = pd.read_csv(value, sep='\t')
                            depth = len(seq_summary.index)
                            Stats(
                                columns=2,
                                items=[
                                    (depth, 'Read count'),
                                ])
                            draw_all_plots(seq_summary, THEME)
        # Insert info table
    if os.stat(args.inserts_json).st_size != 0:
        inserts = json.load(open(args.inserts_json))
        with report.add_section("Insert sequences", "Inserts"):
            seq_segment = pd.read_json(inserts['bed_df'])
            insert_df = seq_segment.merge(stats_df, how='left')
            insert_df['Insert length'] = list(
                map(
                    report_utils.insert_len,
                    insert_df['start'],
                    insert_df['end'],
                    insert_df['Length']))
            insert_df = insert_df.drop(columns='Length')
            p("""
This table shows which primers were found in the consensus sequence
of each sample and where the inserts were found.
""")
            DataTable.from_pandas(insert_df, use_index=False)
        # MSA section if inserts available
        with report.add_section("Multiple Sequence Alignment", "MSA"):
            p("""
This section shows the inserts aligned with each other or a reference
sequence if provided.
""")
            formatted_msa = format_msa(inserts)
            pre(os.linesep.join(formatted_msa))
    if ('OPTIONAL_FILE' not in os.listdir(args.qc_inserts)):
        with report.add_section("Insert variants", "Insert variants"):
            p("""
The following tables and figures are output from bcftools from
finding any variants between
the consensus insert and the provided reference
insert.
""")
            variants_df = report_utils.variant_counts_table(args.qc_inserts, report)
            DataTable.from_pandas(variants_df, use_index=False)
            trans_df = report_utils.trans_counts(args.qc_inserts, report)
            DataTable.from_pandas(trans_df, use_index=False)
    report.write(args.report)
    logger.info(f"Report written to {args.report}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument("report", help="Report output file")
    parser.add_argument(
        "--stats", nargs='*', help="Fastcat per-read stats file(s).")
    parser.add_argument(
        "--versions", required=True,
        help="directory containing CSVs containing name,version.")
    parser.add_argument(
        "--params", default=None, required=True,
        help="A JSON file containing the workflow parameter key/values")
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    parser.add_argument(
        "--per_barcode_stats", nargs='+',
        help="fastcat stats file for each sample before filtering")
    parser.add_argument(
        "--status", nargs='+',
        help="status")
    parser.add_argument(
        "--plannotate_json",
        help="Plannotate Json.")
    parser.add_argument(
        "--inserts_json",
        help="inserts Json.")
    parser.add_argument(
        "--lengths",
        help="report name")
    parser.add_argument(
        "--downsampled_stats",
        nargs='+')
    parser.add_argument(
        "--host_filter_stats", nargs='+',
        help="fastcat stats file after host filtering")
    parser.add_argument(
        "--qc_inserts",
        help="insert vcfs")
    parser.add_argument(
        "--assembly_quality",  nargs='+',
        help="qc quality summaries")
    return parser


if __name__ == "__main__":
    parser = argparser()
    main(parser.parse_args())
