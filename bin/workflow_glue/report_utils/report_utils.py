#!/usr/bin/env python
"""Create tables for the report."""

from dominate.tags import p
from ezcharts.components import bcfstats
from ezcharts.components.ezchart import EZChart
from ezcharts.plots.categorical import barplot
from ezcharts.plots.util import read_files
import pandas as pd

THEME = 'epi2melabs'


def tidyup_status_file(status_sheet, annotations):
    """Tidy up the sample status file."""
    sample_status = pd.read_csv(status_sheet[0], header=None)
    unique_samples = sample_status[0].unique()
    sample_status_dic = {}
    # Default all to success
    for sample in unique_samples:
        sample_status_dic[sample] = 'Completed successfully'
    filter_pass = sample_status[sample_status[1] != 'Completed successfully']
    # Collect failures
    failures = dict(zip(filter_pass[0], filter_pass[1]))
    completed_annotations = list(annotations)
    # If no failure for a sample status sheet then success
    success = sample_status[sample_status[1] == 'Completed successfully']
    # Check for corresponding annotations, if none update status
    no_annotations = success[~success[0].isin(completed_annotations)]
    for sample in list(no_annotations[0]):
        failures[sample] = 'Completed but no annotations found in the database'
    passed_list = unique_samples.tolist()
    # Update sample status dictionary with any failure messages
    # Also create a list of passed samples to iterate through later
    for k, v in failures.items():
        sample_status_dic[k] = v
        if v != 'Completed but failed to reconcile':
            passed_list.remove(k)
    passed_list.sort()
    # Output list of all sample names
    all_sample_names = unique_samples.tolist()
    all_sample_names.sort()
    return (passed_list, all_sample_names, sample_status_dic)


def read_count_barplot(per_barcode_stats, report):
    """Plot per sample read count bar chart."""
    # Per barcode read count
    seq_summary = read_files(per_barcode_stats)
    barcode_counts = (
        pd.DataFrame(seq_summary['sample_name'].value_counts())
        .sort_index()
        .reset_index()
        .rename(
            columns={'index': 'Sample', 'sample_name': 'Count'})
    )
    with report.add_section("Read Counts", "Read Counts"):
        p(
            """Number of reads per sample."""
        )
        order = barcode_counts["Sample"].values.tolist()
        plt = barplot(
            data=barcode_counts, x='Sample', y='Count', order=order)
        plt.xAxis = dict(
            name="Sample",
            type='category',
            axisLabel=dict(interval=0, rotate=45),
            axisTick=dict(alignWithLabel=True))
        plt.yAxis = dict(name="Count", type='value')
        EZChart(
            plt,
            theme='epi2melabs')


def insert_len(start, end, length):
    """Insert length calc."""
    if end >= start:
        insert_length = end - start
    else:
        insert_length = (length - start) + end
    return insert_length


def create_fastcat_dic(sample_names, raw, host_file, downsampled_file):
    """Create dictionary using sample names and fastcat files available."""
    # Collect per sample fastcat stats check no optional files
    if ('host_filter_stats/OPTIONAL_FILE' in host_file):
        host_filt = host_file.remove('host_filter_stats/OPTIONAL_FILE')
        if host_filt is None:
            host_filt = []
    else:
        host_filt = host_file
    if ('downsampled_stats/OPTIONAL_FILE' in downsampled_file):
        summary_stats = downsampled_file.remove(
            'downsampled_stats/OPTIONAL_FILE')
        if summary_stats is None:
            summary_stats = []
    else:
        summary_stats = downsampled_file
    per_sample_dic = {}
    lists = {'Raw': raw, 'Hostfilt': host_filt, 'Downsampled': summary_stats}
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


def variant_counts_table(
        bcf_stats,
        header="**Variant counts:**", report=None):
    """Create a report section contains variant counts.

    :param bcf_stats: one or more outputs from `bcftools stats`.
    :param sample_as_columns: transpose table to put data for each
        sample into a column.
    :param header: a markdown formatted header.
    :param report: an HTMLSection instance.

    :returns: an HTMLSection instance, if `report` was provided the given
        instance is modified and returned.
    """
    bcf_stats = bcfstats.load_bcfstats(bcf_stats)
    p("""
Variant counts per sample. See output bcf
file for info on individual variants.
""")
    df = bcf_stats['SN'].drop(columns='samples')
    return df


def trans_counts(
        bcf_stats,
        header="**Transitions and tranversions:**", report=None):
    """Create a report section with transition and transversion counts.

    :param bcf_stats: one or more outputs from `bcftools stats`.
    :param header: a markdown formatted header.
    :param report: an HTMLSection instance.

    :returns: an HTMLSection instance, if `report` was provided the given
        instance is modified and returned.
    """
    bcf_stats = bcfstats.load_bcfstats(bcf_stats)
    p("""
Trans counts per sample.
See output bcf file for info on individual transitions.
""")
    df = bcf_stats['TSTV']
    return df
