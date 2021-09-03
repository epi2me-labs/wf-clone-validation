#!/usr/bin/env python

"""Report component for displaying information from wf-clone-validation."""
import argparse
import json
import os


from aplanat import bars, report
from aplanat import json_item
from aplanat.components import fastcat
from aplanat.components import simple as scomponents
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
from spoa import poa


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


def per_assembly(assemblies_file, database, sample_files, item):
    """Build_per_assembly_tuples."""
    plot, annotations = run_plannotate(sample_files['fasta'], database)
    df = pd.read_csv(assemblies_file, sep='\t', index_col=0)
    filename = str(item) + '.final.fasta'
    avg_len = df.loc[filename, 'avg_len']
    tup = {'sample_name': item,
           'plot': plot,
           'annotations': annotations,
           'avg_len': avg_len}
    return(tup)


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
    status_df = pd.DataFrame(pass_fail_dic.items(),
                             columns=['Sample', 'pass/failed reason'])
    sort_df = status_df['Sample'].astype(str).argsort()
    status_df = status_df.iloc[sort_df]
    status_df.to_csv('sample_status.txt', index=False)
    return(status_df, passed_list, all_sample_names, pass_fail_dic)


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


def reverse_complement(seq):
    """Read a seq return reverse complement."""
    comp = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'X': 'X', 'N': 'N'}
    comp_trans = seq.maketrans(''.join(comp.keys()), ''.join(comp.values()))
    return seq.translate(comp_trans)[::-1]


def read_seqkit(bed_files, sep='\t'):
    """Read seqkit bed files and join to single dataframe."""
    dfs = list()
    bed_dic = {}
    for fname in sorted(bed_files):
        df = pd.read_csv(fname, sep=sep,
                         header=None, names=(
                             ['Sample', 'start',
                              'end', 'primer', 'score',
                              'strand', 'sequence']))
        seq = str(df['sequence'][0])
        # need to find the sequence if seq spans origin
        if seq == 'nan':
            file_name = os.path.join(
                        'assemblies/', df['Sample'][0] + '.final.fasta')
            with open(file_name, "r") as fp:
                whole_seq = fp.readlines()[1][:-1:]
            rev_comp = reverse_complement(whole_seq)
            strand_seq = {'-': rev_comp, '+': whole_seq}
            parse_seq = strand_seq[str(df['strand'][0])]
            final_seq = parse_seq[df['start'][0]::] + parse_seq[:df['end'][0]:]
            df['sequence'][0] = final_seq
        else:
            pass
        bed_dic[df['Sample'][0]] = df['sequence'][0]
        dfs.append(df)
    bed_df = pd.concat(dfs).drop(['score', 'sequence'], axis=1)
    return (bed_df, bed_dic)


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
    exec_summary.append('No. reads',
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
    parser.add_argument(
        "--primer_beds", nargs='+',
        help="bed files of extracted sequences")
    parser.add_argument(
        "--align_ref", nargs='+',
        help="insert alignment reference file")
    parser.add_argument(
        "--params", default=None,
        help="A csv containing the parameter key/values")
    parser.add_argument(
        "--versions",
        help="directory contained CSVs containing name,version.")
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
    placeholder = report_doc.add_section(key='stats_table')
    section = report_doc.add_section()
    section.markdown("""
## Assemblies
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
    sample_stats = []
    plannotate_collection = {}
    for item in sample_names:
        if item in passed_samples:
            sample_files = sample_data[item]
            tup_dic = per_assembly(assembly, database, sample_files, item)
            final_samples.append(tup_dic)
            section = report_doc.add_section()
            section.markdown('## Sample: {}'.format(str(item)))
            section.markdown('####{}'.format(pass_fail[3][item]))
            fast_cat_tabs = fast_cat_dic[item]
            alltabs = []
            for key, value in fast_cat_tabs.items():
                alltabs.append(fastcat_report_tab(value, key))
            cover_panel = Tabs(tabs=alltabs)
            plasmidplot = [tup_dic['plot']]
            stats_table = [item] + [int(tup_dic['avg_len'])]
            sample_stats.append(stats_table)
            section.plot(layout(
                            [[cover_panel, plasmidplot]],
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
            section.markdown('## Sample Failed: {}'.format(str(item)))
            section.markdown('#{}'.format(pass_fail[3][item]))
            fast_cat_tabs = fast_cat_dic[item]
            alltabs = []
            for key, value in fast_cat_tabs.items():
                alltabs.append(fastcat_report_tab(value, key))
            cover_panel = Tabs(tabs=alltabs)
            section.plot(cover_panel)
            stats_table = [item] + ['N/A']
            sample_stats.append(stats_table)
    json_object = json.dumps(plannotate_collection, indent=4)
    json_file.write(json_object)
    json_file.close()
    output_feature_table(final_samples)
    section = report_doc.add_section()
    status = args.status
    pass_fail = tidyup_status_file(status)
    placeholder.markdown("##Sample status")
    stats_df = pd.DataFrame(sample_stats, columns=["Sample",
                                                   "Length"])
    merged_inner = pd.merge(pass_fail[0], stats_df)
    placeholder.table(merged_inner, index=False, key='stats_table')
    section = report_doc.add_section()
    # find inserts and put in dir
    current_directory = os.getcwd()
    make_directory = os.path.join(current_directory, r'inserts')
    if not os.path.exists(make_directory):
        os.makedirs(make_directory)
    if args.primer_beds[0] != "primer_beds/OPTIONAL_FILE":
        seq_segment = read_seqkit(args.primer_beds)[0]
        section.markdown("""
### Insert sequences
This table shows which primers were found in the consensus sequence
of each sample and where the inserts were found.
""")
        section.table(seq_segment, key='sequences')
        inserts = ''
        allseq = []
        names = []
        # Align with reference if included
        if args.align_ref[0] != 'OPTIONAL_FILE':
            with open(args.align_ref[0]) as f:
                ref_seq = f.readline().strip()
                allseq.append(ref_seq)
                names.append('Reference')
        for k, v in read_seqkit(args.primer_beds)[1].items():
            inserts += str(k) + ' ' + str(v) + '</br>'
            insert_fn = os.path.join('inserts/', str(k) + '.insert.fasta')
            with open(insert_fn, "a") as fp:
                fp.write('>' + k + '\n' + v + '\n')
            allseq.append(v)
            names.append(k)
        msa_report = []
        msa = poa(allseq)[1]
        section.markdown("""
### Multiple sequence alignment
This section shows the inserts aligned with each other or a reference
sequence if provided.
""")
        # make sure names are all same length for MSA
        msa_names = []
        for name in names:
            new_name = name + (' '*(max(map(len, names))-len(name)))
            msa_names.append(new_name)
        for i in range(0, len(msa)):
            msa_report.append(msa_names[i] + ' ' + msa[i])
        section.markdown("<pre>" + os.linesep.join(msa_report) + "</pre>")

    # Versions and params
    report_doc.add_section(
        section=scomponents.version_table(args.versions))
    report_doc.add_section(
        section=scomponents.params_table(args.params))
    report_doc.write('wf-clone-validation-report.html')


if __name__ == "__main__":
    main()
