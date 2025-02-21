#!/usr/bin/env python
"""Create tables for the report."""

from bokeh.models import ColumnDataSource, Segment
from bokeh.plotting import figure
from dominate.tags import p
from ezcharts.components import bcfstats
import pandas as pd


THEME = 'epi2melabs'


def tidyup_status_file(status_sheet, annotations):
    """Tidy up the sample status file."""
    sample_status = pd.read_csv(status_sheet[0], header=None, dtype={0: str})
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
    lists = {'Raw': raw, 'Host_filtered': host_filt, 'Downsampled': summary_stats}
    for sample in sample_names:
        new_dic = {}
        item_search = '/' + str(sample) + '.'
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


def dotplot_assembly(maf_assembly):
    """Dotplot of assembly using a .maf format file."""
    # Create a bokeh plot
    plot = figure(
        width=400, height=400,
        min_border=0, x_axis_label="position",
        y_axis_label="position", title="Dot plot")
    plot.toolbar_location = None
    # Iterate through maf file to create a dataframe
    records = list()
    with open(f"mafs/{maf_assembly}") as maf:
        while True:
            line = maf.readline()
            if line.startswith('#'):
                # get read length
                if 'letters' in line:
                    read_length = int(line.split('letters=')[1])
                    continue
                else:
                    continue
            # a is for each alignment
            elif line.startswith('a'):
                # take successive 's' lines
                r1 = maf.readline().split()[1:5]
                r2 = maf.readline().split()[1:5]
                maf.readline()
                records.append(r1 + r2)
            elif line == "":
                break
            else:
                raise IOError("Cannot read alignment file")
    # take reference start, length and orientation
    # and query start, length and orientation
    names = [
        'ref', 'rstart', 'rlen', 'rorient',
        'query', 'qstart', 'qlen', 'qorient']
    df_all = pd.DataFrame(records, columns=names)
    df_all = df_all.astype({'qstart': int, 'qlen': int, 'rstart': int, 'rlen': int})
    # If query orientation is +
    df = df_all[df_all.qorient.isin(['+'])]
    # create query and ref end columns by adding length to start
    df['qend'] = df['qstart'] + df['qlen']
    df.loc[df['rorient'] == '+', 'rend'] = df['rstart'] + df['rlen']
    # If reference orientation is negative switch reference start and end
    df.loc[df['rorient'] == '-', 'rend'] = df['rstart']
    df.loc[df['rorient'] == '-', 'rstart'] = df['rstart'] - df['rlen']
    # Add fwd lines to plot
    source = ColumnDataSource(df)
    glyph = Segment(x0='rstart', y0='qstart', x1='rend', y1='qend', line_color="black")
    plot.add_glyph(source, glyph)
    # if query orientation is -
    df = df_all[df_all.qorient.isin(['-'])]
    # If the orientation is "-", start coordinate is in the reverse strand (maf docs)
    # Therefore as plot will be + vs + query start needs to be flipped
    df['qstart'] = read_length - df['qstart']
    df['qend'] = df['qstart'] - df['qlen']
    df['rend'] = df['rstart'] + df['rlen']
    # Add reverse complement lines to plot
    source = ColumnDataSource(df)
    glyph = Segment(x0='rstart', y0='qstart', x1='rend', y1='qend', line_color="red")
    plot.add_glyph(source, glyph)
    return plot


def bamstats_table(input_files, passed_samples):
    """Use bamstats files to create table with ref and identity."""
    df = pd.concat(
        (pd.read_csv(f, sep='\t') for f in input_files), ignore_index=True)
    df = df[['sample_name', 'ref_coverage', 'coverage', 'acc']]
    df = df.rename(
        columns={
            "sample_name": "Sample name",
            "ref_coverage": "Reference coverage",
            "coverage": "Assembly coverage",
            "acc": "BLAST Identity"}
    )
    # Add 0's for any samples which have an assembly but no alignments.
    df_passed = pd.DataFrame(passed_samples, columns=['Sample name'])
    df = pd.merge(df_passed, df, on='Sample name', how='outer').fillna(0)
    return df


def get_cutsite_table(cutsite_csv, sample_names):
    """Use cutsite csv to create linearisation efficiency table."""
    # Make a default dataframe with all samples
    default_df = pd.DataFrame(sample_names, columns=['Sample'])
    default_df['Linearisation efficiency (%)'] = "N/A"
    # Create a dataframe from the cutsite_csv
    cutsite_table = pd.read_csv(
        cutsite_csv,
        index_col=False,
        names=['Sample', 'readcount', 'cutsitecount'],
        header=None)
    cutsite_table['Linearisation efficiency (%)'] = (
        (cutsite_table['readcount'] -
            cutsite_table['cutsitecount']) /
        cutsite_table['readcount']) * 100
    cutsite_table = cutsite_table.round(2)
    cutsite_table = cutsite_table.drop(
        ['readcount', 'cutsitecount'], axis=1)
    # Update the default table where values available
    default_df.update(cutsite_table)
    cutsite_table = default_df.sort_values(by='Sample')
    return cutsite_table
