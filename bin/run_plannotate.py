#!/usr/bin/env python

"""Report component for displaying information from wf-clone-validation."""
import argparse
import json
import os


from aplanat import json_item
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
    old_df = df.copy(deep=True)
    clean_df = clean_results(df)
    old_df.reset_index(drop=True, inplace=True)
    print(clean_df, old_df)
    return plot, clean_df, old_df


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


def bed_file(item, df):
    """Bed format for annotations."""
    display_columns = [
        'Start Location', 'End Location',
        'Feature',
        'Strand']
    df = df[display_columns]
    df["Strand"] = df["Strand"].apply(pd.to_numeric)
    df.loc[df['Strand'] == 0, 'Strand'] = "-"
    df.loc[df['Strand'] == 1, 'Strand'] = "+"
    df.insert(0, 'Name', value=str(item))
    df.to_csv(
        str(item)+'.annotations.bed', sep="\t", header=False, index=False)
    return df


def per_assembly(database, sample_file, item):
    """Run plannotate for a sample.

    :param database:
    :param sample_file:
    :param item: the sample
    """
    plot, annotations, clean_df = run_plannotate(sample_file, database)
    bed_file(item, annotations)
    with pysam.FastxFile(sample_file) as fh:
        seq_len = len(next(fh).sequence)
    tup = {
        'sample_name': item,
        'plot': plot,
        'annotations': annotations,
        'seq_len': seq_len}
    return(tup, clean_df)


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


def main():
    """Entry point to create a wf-clone-validation report."""
    parser = argparse.ArgumentParser(
        'Plannotate annotation',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False)
    parser.add_argument(
        "--database", default='unknown',
        help="database to use, directory containing BLAST et. al. files.")
    parser.add_argument(
        "--sequences",
        help="sequences in directory to run plannotate on.",
        required=False)
    args = parser.parse_args()

    final_samples = []
    report_dic = {}
    plannotate_collection = {}
    json_file = open("plannotate.json", "a")
    if args.sequences:
        for filename in os.listdir(args.sequences):
            sample_file = os.path.join(args.sequences, filename)
            name = str(filename).split('.')[0]
            tup_dic, clean_df = per_assembly(args.database, sample_file, name)
            final_samples.append(tup_dic)
            plasmid_len = tup_dic['annotations']['Plasmid length'][0]
            feature_dic = tup_dic['annotations'].drop(
                            ['Plasmid length'], axis=1)
            features = feature_dic.to_dict('records')
            output_json = json_item(tup_dic['plot'])
            plannotate_dic = {
                "reflen": float(plasmid_len),
                "features": features,
                "plot": output_json['doc']}
            plannotate_collection[name] = plannotate_dic
            report = {}
            report['sample_name'] = name
            report['plot'] = clean_df.to_json()
            report['annotations'] = tup_dic['annotations'].to_json()
            report['seq_len'] = tup_dic['seq_len']
            report_dic[name] = report

    # outputs for epi2me
    output_feature_table(final_samples)
    json_object = json.dumps(plannotate_collection, indent=4)
    json_file.write(json_object)
    json_file.close()

    # outputs for report
    json_file = open("plannotate_report.json", "a")
    json_object = json.dumps(report_dic)
    json_file.write(json_object)
    json_file.close()


if __name__ == "__main__":
    main()
