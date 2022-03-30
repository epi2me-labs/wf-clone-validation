#!/usr/bin/env python

"""Finding inserts."""

import argparse
import json
import os

import pandas as pd
from spoa import poa


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
        df = pd.read_csv(
                fname, sep=sep,
                header=None, names=(
                    [
                        'Sample', 'start',
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
    return bed_df, bed_dic


def make_msa(inserts_dic, *reference):
    """Make multiple sequence alignment."""
    allseq = []
    names = []
    # Align with reference if included
    if reference:
        with open(reference[0]) as f:
            ref_seq = f.readline().strip()
            allseq.append(ref_seq)
            names.append('Reference')
    for k, v in inserts_dic.items():
        allseq.append(v)
        names.append(k)
    msa_report = []
    msa = poa(allseq)[1]
    # make sure names are all same length for MSA
    msa_names = []
    for name in names:
        new_name = name + (' '*(max(map(len, names))-len(name)))
        msa_names.append(new_name)
    for i in range(0, len(msa)):
        msa_report.append(msa_names[i] + ' ' + msa[i])
    return(msa_report)


def main():
    """Entry point to create a wf-clone-validation report."""
    parser = argparse.ArgumentParser(
        'Clone Validation QC report',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False)
    parser.add_argument(
        "--primer_beds", nargs='+',
        help="bed files of extracted sequences",
        required=False)
    parser.add_argument(
        "--reference", nargs='+',
        help="reference", required=False)
    args = parser.parse_args()
    f = open('insert_data.json', 'w')
    if args.primer_beds:
        # find inserts and put in dir
        current_directory = os.getcwd()
        make_directory = os.path.join(current_directory, r'inserts')
        if not os.path.exists(make_directory):
            os.makedirs(make_directory)
        inserts = ''
        read_seqkit(args.primer_beds)
        seqkit = read_seqkit(args.primer_beds)
        for k, v in seqkit[1].items():
            inserts += str(k) + ' ' + str(v) + '</br>'
            insert_fn = os.path.join('inserts/', str(k) + '.insert.fasta')
            with open(insert_fn, "a") as fp:
                fp.write('>' + k + '\n' + v + '\n')
        if args.reference:
            msa = make_msa(seqkit[1], args.reference)
        else:
            msa = make_msa(seqkit[1])
        inserts_json = {
            'bed_df': seqkit[0].to_json(), 'bed_dic': seqkit[1], 'msa': msa}
        json.dump(inserts_json, f)


if __name__ == "__main__":
    main()
