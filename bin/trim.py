#!/usr/bin/env python
"""Trim Canu assembly sequences."""

import argparse
import sys

import pysam


def get_output_handler(path):
    """Open file path or stdout."""
    if path != '-':
        fh = open(path, 'w')
    else:
        fh = sys.stdout
    return fh


def trim(entry):
    """Trim following fastq comment suggestion."""
    # example Canu output fastq header >tig00000001 len=3207 reads=43 class=contig
    # suggestRepeat=no suggestBubble=no suggestCircular=no trim=0-3207
    split = {}
    for i in entry.comment.split(' '):
        subsplit = i.split('=')
        split[subsplit[0]] = subsplit[1]

    trim = [int(i) for i in split['trim'].split('-')]
    trimmed = entry.sequence[trim[0]:trim[1]]

    return trimmed


def main(args):
    """For each sequence, trim and write to output."""
    handler = get_output_handler(args.output)
    for entry in pysam.FastxFile(args.sequence):
        if 'trim=' not in entry.comment:
            continue
        handler.write(f">{entry.name}\n{trim(entry)}\n")
    handler.close()


def argparser():
    """Argument parser for entrypoint."""
    parser = argparse.ArgumentParser(
        "trim",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False)
    parser.add_argument(
        dest="sequence",
        help="File in .FASTA format containing a single sequence/contig."
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        default="-",
        help="Path at which to write the fixedsequence/contig.",
        required=False
    )
    return parser


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
