#!/usr/bin/env python

"""Finding inserts."""

import argparse
import json
import os

import pandas as pd
from pysam import FastaFile
from spoa import poa
from .util import wf_parser  # noqa: ABS101


def reverse_complement(seq):
    """Read a seq return reverse complement."""
    comp = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'X': 'X', 'N': 'N'}
    comp_trans = seq.maketrans(''.join(comp.keys()), ''.join(comp.values()))
    return seq.translate(comp_trans)[::-1]


def read_seqkit(bed_files, assemblies_dir):
    """Read beds and return single data frame and dict of the insert seqs."""
    dfs = list()
    bed_dic = {}
    for fname in sorted(bed_files):
        df = pd.read_csv(
                fname, sep='\t',
                header=None, names=(
                    [
                        'Sample', 'start',
                        'end', 'primer', 'score',
                        'strand', 'sequence']),
                dtype={'Sample': str})
        seq = str(df['sequence'][0])
        # need to find the sequence if seq spans origin
        if seq == 'nan':
            file_name = os.path.join(
                assemblies_dir, df['Sample'][0] + '.final.fasta')
            with open(file_name, "r") as fp:
                whole_seq = fp.readlines()[1][:-1:]
            final_seq = whole_seq[df['start'][0]::] + whole_seq[:df['end'][0]:]
            if df['strand'][0] == '-':
                final_seq = reverse_complement(final_seq)
            df['sequence'][0] = final_seq
        else:
            pass
        bed_dic[df['Sample'][0]] = df['sequence'][0]
        dfs.append(df)
    bed_df = pd.concat(dfs).drop(['score', 'sequence'], axis=1)
    bed_df.reset_index(drop=True, inplace=True)
    return bed_df, bed_dic


def make_msa(inserts_dic, reference=None):
    """Make multiple sequence alignment."""
    allseq = []
    names = []
    # Align with reference if included
    if reference:
        refseq = FastaFile(reference)
        ref_seq = refseq.fetch(refseq.references[0])
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
    return (msa_report)


def main(args):
    """Entry point to create a wf-clone-validation report."""
    with open(args.output, 'w') as f:
        ref = "No reference"
        if args.reference:
            insert_ref = FastaFile(args.reference).references
            if len(insert_ref) > 1:
                raise ValueError(
                    f"""Insert reference can only contain one fasta record;
                     {len(insert_ref)} found.""")
            ref = insert_ref[0]
        if args.primer_beds:
            # find inserts and put in dir
            current_directory = os.getcwd()
            make_directory = os.path.join(current_directory, r'inserts')
            # only create inserts directory if there are primer_beds
            # as optional output in nextflow
            if not os.path.exists(make_directory):
                os.makedirs(make_directory)
            sk_df, sk_dict = read_seqkit(args.primer_beds, args.assemblies)
            for k, v in sk_dict.items():
                insert_fn = os.path.join(args.insert_dir, f'{str(k)}.insert.fasta')
                with open(insert_fn, "a") as fp:
                    fp.write('>' + str(k) + '\n' + str(v) + '\n')
            # If assembly will be large, skip MSA creation
            if args.large_construct:
                inserts_json = {
                    'bed_df': sk_df.to_json(),
                    'bed_dic': sk_dict,
                    'reference': ref}
            else:
                # If reference is available, it will be included to perform the
                # Multiple sequence alignment.
                # Otherwise MSA will be done using the inserts available
                if args.reference:
                    msa = make_msa(sk_dict, args.reference)
                else:
                    msa = make_msa(sk_dict)
                inserts_json = {
                    'bed_df': sk_df.to_json(),
                    'bed_dic': sk_dict, 'msa': msa,
                    'reference': ref}
            json.dump(inserts_json, f)


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("find_inserts")
    parser.add_argument(
        "--output",
        help="output json name",
        required=True
    )
    parser.add_argument(
        "--primer_beds", nargs='+',
        help="bed files of extracted sequences",
        required=False)
    parser.add_argument(
        "--large_construct", default=False, action="store_true",
        help="large construct mode skip msa",
        required=False)
    parser.add_argument(
        "--reference",
        help="reference", required=False)
    parser.add_argument(
        "--insert_dir", default="inserts",
        help="output directory for insert fastas"
    )
    parser.add_argument(
        "--assemblies", default="assemblies",
        help="Full assemblies directory for finding insert"
    )
    return parser


if __name__ == "__main__":
    args = argparse().parse_args()
    main(args)
