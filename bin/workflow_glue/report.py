#!/usr/bin/env python
"""Create workflow report."""
from itertools import count, takewhile
import json
import os

from dominate.tags import p, pre
from dominate.util import raw
from ezcharts.components.ezchart import EZChart
from ezcharts.components.fastcat import load_histogram, load_stats, SeqSummary
from ezcharts.components.plotmetadata import read_count_barplot
from ezcharts.components.reports import labs
from ezcharts.layout.snippets import Stats, Tabs
from ezcharts.layout.snippets.table import DataTable
from ezcharts.plots import BokehPlot
from ezcharts.plots.util import read_files
from ezcharts.util import get_named_logger  # noqa: ABS101
import numpy as np
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


def format_badge(status):
    """Format badge, green for success, red for failure."""
    if status == "Completed successfully" or status is True:
        badge = 'badge-pass-solid'
        if status is True:
            status = 'Pass'
    else:
        badge = 'badge-fail-solid'
        if not isinstance(status, str):
            status = 'Fail'
    return '<span class="badge badge-icon-solid rounded-pill p-2 ' + \
        badge + '">' + status + '</span>'


def add_expected_column(
        bam_stats_df, sample_status_df, column_name,
        expected_coverage, expected_identity):
    """Add expected column to the sample status using BAM stats."""
    bam_stats_df = bam_stats_df.copy()
    bam_stats_df["expected"] = np.where(
        (bam_stats_df['Reference coverage'] >= expected_coverage) &
        (bam_stats_df['Assembly coverage'] >= expected_coverage) &
        (bam_stats_df['BLAST Identity'] >= expected_identity),
        True, False)
    merged = pd.merge(sample_status_df, bam_stats_df[
                ['Sample name', 'expected']],
                how="outer", left_on='Sample', right_on='Sample name')
    merged[column_name] = merged.apply(
            lambda x: format_badge(x['expected']),
            axis=1)
    merged = merged.drop('Sample name', axis=1)
    merged = merged.drop('expected', axis=1)
    return merged


def main(args):
    """Run the entry point."""
    logger = get_named_logger("Report")
    report = labs.LabsReport(
        "Clone validation report", "wf-clone-validation-new",
        args.params, args.versions, args.wf_version)
    # Define criteria for expected assembly tick in status table
    expected_coverage = args.expected_coverage
    expected_identity = args.expected_identity
    plannotate_annotations = json.load(open(args.plannotate_json)).keys()
    # Get sample info from status file
    passed_samples, sample_names, sample_status_dic = \
        report_utils.tidyup_status_file(
            args.status, plannotate_annotations)
    # Sample status table place holder
    sample_status_table = report.add_section("Sample status", "Sample status")
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

The plannotate plot may have overlapping annotation labels,
use the zoom and hover tools to decipher the labels.
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
                    EZChart(
                        bk_plot, 'epi2melabs',
                        width='100%', height='100%')
                    DataTable.from_pandas(annotations, use_index=False)
    # Per barcode read count plot
    with report.add_section("Read Counts", "Read Counts"):
        p(
            """Number of reads per sample."""
        )
        read_count_barplot(args.metadata)
    with report.add_section("Read stats", "Read stats"):
        p("""
For each assembly, read length statistics and plots of quality \
before and after host filtering, and after downsampling \
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
                        with internal_tabs.add_tab(key.replace("_", " ")):
                            if key == "Raw":
                                stats = load_histogram(value)
                                depth = stats["count"].sum()
                            else:
                                stats = load_stats(value)
                                depth = len(stats.index)
                            Stats(
                                columns=2,
                                items=[
                                    (f"{depth:,d}", 'Read count'),
                                ])
                            SeqSummary(value)
    # Insert info table
    json_combined = {}
    for inserts_json in args.inserts_json:
        # the inserts JSON files are empty if the wf was run without primers
        if os.stat(inserts_json).st_size != 0:
            inserts_data = json.load(open(inserts_json))
            json_combined[inserts_data["reference"]] = inserts_data
    if json_combined:
        with report.add_section("Insert sequences", "Inserts"):
            tabs = Tabs()
            with tabs.add_dropdown_menu():
                for ref, inserts in json_combined.items():
                    with tabs.add_dropdown_tab(ref):
                        seq_segment = pd.read_json(
                            inserts['bed_df'], dtype={'Sample': str})
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
        # MSA section if available
        msa_jsons = {r: i for r, i in json_combined.items() if "msa" in i}
        if msa_jsons:
            with report.add_section("Multiple Sequence Alignment", "MSA"):
                p("""
                This section shows the inserts aligned with each other or a reference
                sequence if provided.
                """)
                tabs = Tabs()
                with tabs.add_dropdown_menu():
                    for ref, inserts in msa_jsons.items():
                        with tabs.add_dropdown_tab(ref):
                            formatted_msa = format_msa(inserts)
                            pre(os.linesep.join(formatted_msa))

    # Dictionary to collect fields to be added to sample status table
    sample_status_fields = {}
    if ('OPTIONAL_FILE' not in os.listdir(args.qc_inserts)):
        if args.insert_alignment_bamstats:
            insert_qc_df = report_utils.bamstats_table(
                args.insert_alignment_bamstats,
                passed_samples)
            sample_status_fields["Expected insert"] = insert_qc_df
            with report.add_section("Insert QC", "Insert QC"):
                p("""
This section can be used to ensure the provided insert reference
matches the insert found in the assembly.

The table below shows coverage and BLAST identity between
the two.

Reference coverage is the percentage of the provided insert reference
sequence covered in the alignment with the assembled construct.

Assembly coverage is the percentage of the assembled insert
sequence covered in the alignment with the provided insert reference.

BLAST identity is calculated as: (length - ins - del - sub) / length.

If both coverage and identity are 0,
the assembled insert did not align with the provided insert
reference.
""")
                DataTable.from_pandas(insert_qc_df, use_index=False)

        with report.add_section("Insert variants", "Insert variants"):
            p("""
        The following tables and figures are output from bcftools from
        finding any variants between
        the consensus insert and the provided reference
        insert.
            """)
            variants_df = report_utils.variant_counts_table(args.qc_inserts, report)
            variants_df = variants_df.sort_values(by='id')
            DataTable.from_pandas(variants_df, use_index=False)
            trans_df = report_utils.trans_counts(args.qc_inserts, report)
            trans_df = trans_df.sort_values(by='id')
            DataTable.from_pandas(trans_df, use_index=False)
    # Full plasmid QC section
    # Handling for if no assemblies have aligned to the reference.
    full_ref_analysis = False
    with open(args.metadata) as data_file:
        data = json.load(data_file)
        for item in data:
            if "full_reference" in item:
                full_ref_analysis = True
                break
    if passed_samples and full_ref_analysis:
        if not args.reference_alignment_bamstats:
            with report.add_section("Full construct QC", "Construct QC"):
                p("""
None of the assemblies found aligned with the provided reference.
""")
    if args.reference_alignment_bamstats:
        # Use bamstats to output table with coverage and identity per sample
        assembly_df = report_utils.bamstats_table(
            args.reference_alignment_bamstats,
            passed_samples)
        sample_status_fields["Expected Assembly"] = assembly_df
        with report.add_section("Full construct QC", "Construct QC"):
            p("""
This section can be used to ensure the provided reference matches the assembly.

The table belows shows coverage and BLAST identity between
the provided reference and assembly.

Reference coverage is the percentage of the provided reference sequence covered
in the alignment with the assembled construct.

Assembly coverage is the percentage of assembled
construct sequence covered in the alignment with the provided reference.

BLAST Identity is calculated as: (length - ins - del - sub) / length.

If both coverage and identity are 0, the assembly did not align with the provided
reference.
""")
            DataTable.from_pandas(assembly_df, use_index=False)
            if args.full_assembly_variants:
                p("""
    Additionally, BCFtools was used to report any variants between
    the provided reference and assembly.
    """)
                # Use BCF stats report to output table
                # with summary of per sample variant counts.
                variants_df = report_utils.variant_counts_table(
                    args.full_assembly_variants, report)
                variants_df = variants_df.sort_values(by='id')
                DataTable.from_pandas(variants_df, use_index=False)
                # Also use VCF stats report to output table
                # with summary of per sample transition and transversion counts.
                trans_df = report_utils.trans_counts(
                    args.full_assembly_variants, report)
                trans_df = trans_df.sort_values(by='id')
                DataTable.from_pandas(trans_df, use_index=False)
    # Add info to the status table in the sample status section
    with sample_status_table:
        p("""
This table gives the status of each sample, either completed \
successfully or the reason the assembly failed. \
If applicable the mean quality for the whole \
construct has been provided, derived by Medaka \
from aligning the reads to the consensus.
""")
        if args.assembly_tool == "flye":
            raw("""
The assembly was generated using <a href="https://github.com/fenderglass/Flye">Flye</a>.
        """)
        else:
            # Must be canu
            raw("""
The assembly was generated using <a href="https://github.com/marbl/canu">Canu</a>.
        """)
        if args.reference_alignment_bamstats or args.insert_alignment_bamstats:
            # If references provided add explanation of the expected_assembly
            # and expected insert columns to the report.
            raw("""
The assemblies and/or inserts were aligned with provided references and
marked as expected if they meet both the acceptance
criteria defined by the expected_coverage and expected_identity
parameters which have been set to {}% and {}% respectively.
""".format(expected_coverage, expected_identity))

        status_df = pd.DataFrame(
            sample_status_dic.items(),
            columns=['Sample', 'Assembly completed / failed reason'])
        sort_df = status_df['Sample'].astype(str).argsort()
        status_df = status_df.iloc[sort_df]
        # Length
        merged_status_df = pd.merge(status_df, stats_df)
        merged_status_df.to_csv('sample_status.txt', index=False)
        # Mean quality
        if ('assembly_quality/OPTIONAL_FILE' not in args.assembly_quality):
            qc_df = read_files(args.assembly_quality, dtype={'sample_name': str})[
                ['sample_name', 'mean_quality']]
            qc_df = qc_df.rename(columns={
                'sample_name': 'Sample', 'mean_quality': 'Mean Quality'})
            qc_df = qc_df.reset_index(drop=True)
            merged_status_df = pd.merge(merged_status_df, qc_df, how="outer")
            # N/A for samples which failed assembly
            merged_status_df.fillna('N/A', inplace=True)
        merged_status_df['Assembly completed / failed reason'] = \
            merged_status_df.apply(
                lambda x: format_badge(
                    x['Assembly completed / failed reason']),
                axis=1)
        # Expected assembly/inserts
        for expected_column in sample_status_fields.keys():
            merged_status_df = add_expected_column(
                sample_status_fields[expected_column], merged_status_df,
                expected_column,
                expected_coverage, expected_identity)
        DataTable.from_pandas(merged_status_df, use_index=False)
    # dot plots
    if ('OPTIONAL_FILE' not in os.listdir(args.mafs)):
        with report.add_section("Dot plots", "Dot plots"):
            raw("""
These dot plots have been created by aligning the assembly to itself
to reveal any repeats or repetitive regions.
Black is for repeats found in the forward strand and
red for repeats found in the reverse complement.
This was done using
<a href=
"https://github.com/UCSantaCruzComputationalGenomicsLab/last/tree/master">
last</a>
""")
            tabs = Tabs()
            with tabs.add_dropdown_menu():
                for item in sorted(os.listdir(args.mafs)):
                    sample_name = item.split('.')[0]
                    with tabs.add_dropdown_tab(str(sample_name)):
                        bk_plot = BokehPlot()
                        bk_plot._fig = report_utils.dotplot_assembly(item)
                        EZChart(
                            bk_plot, 'epi2melabs', height="500px", width="500px")

    client_fields = None
    if args.client_fields:
        with open(args.client_fields) as f:
            try:
                client_fields = json.load(f)
            except json.decoder.JSONDecodeError:
                error = "ERROR: Client info is not correctly formatted"

        with report.add_section("Workflow Metadata", "Workflow Metadata"):
            raw("""This table records the data specified in the client fields
                input.
                """)

            if client_fields:
                df = pd.DataFrame.from_dict(
                    client_fields, orient="index", columns=["Value"])
                df.index.name = "Key"

                # Examples from the client had lists as values so join lists
                # for better display
                df['Value'] = df.Value.apply(
                    lambda x: ', '.join(
                        [str(i) for i in x]) if isinstance(x, list) else x)

                DataTable.from_pandas(df)
            else:
                p(error)
    if args.cutsite_csv:
        with report.add_section("Linearisation efficiency", "Linearisation"):
            p("""This table gives the percentage of reads which did not align
                across the cutsite provided in the sample sheet, and are therefore
                assumed to be linearised correctly.\n""")
            p("""This is calculated by comparing the reference sequence
                with the reads, so this metric will still be
                produced for failed assemblies.\n""")
            p("""Furthermore, if a cut site was not provided for a
                sample, this metric will be N/A.\n""")
            cutsite_table = report_utils.get_cutsite_table(
                args.cutsite_csv, sample_names)
            DataTable.from_pandas(cutsite_table, use_index=False)
    report.write(args.report)
    logger.info(f"Report written to {args.report}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument("report", help="Report output file")
    parser.add_argument(
        "--metadata", default='metadata.json', required=True,
        help="sample metadata")
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
        "--inserts_json", nargs='+',
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
    parser.add_argument(
        "--mafs",
        help="mafs")
    parser.add_argument(
        "--cutsite_csv",
        help="cutsite csv file")
    parser.add_argument(
        "--full_assembly_variants",
        help="Enable BCF stats reports for full-plasmid reference")
    parser.add_argument(
        "--assembly_tool", choices=['canu', 'flye'],
        help="Assembly tool selected Canu or Flye")
    parser.add_argument(
        "--reference_alignment_bamstats", nargs='+',
        help="bamstats files from reference alignment")
    parser.add_argument(
        "--insert_alignment_bamstats", nargs='+',
        help="bamstats files from insert alignment")
    parser.add_argument(
        "--client_fields",
        help="JSON file containing client_fields")
    parser.add_argument(
        "--wf_version", default='unknown',
        help="version of the executed workflow")
    parser.add_argument(
        "--expected_coverage", required=True,
        type=float,
        help="Expected coverage as a percentage")
    parser.add_argument(
        "--expected_identity", required=True,
        type=float,
        help="Expected identity as a percentage")

    return parser


if __name__ == "__main__":
    parser = argparser()
    main(parser.parse_args())
