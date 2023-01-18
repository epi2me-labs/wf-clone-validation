#!/usr/bin/env python
"""Go deconcatenate your sequences."""

import argparse
import sys

import pysam
from .util import wf_parser  # noqa: ABS101


def get_output_handler(path):
    """Open file path or stdout."""
    if path != '-':
        fh = open(path, 'w')
    else:
        fh = sys.stdout
    return fh


def trim(entry):
    """Trim following fastq comment suggestion."""
    split = {}
    for i in entry.comment.split(' '):
        subsplit = i.split('=')
        split[subsplit[0]] = subsplit[1]

    trim = [int(i) for i in split['trim'].split('-')]

    if len(trim) > 1:
        trimmed = entry.sequence[trim[0]:trim[1]]
    else:
        trimmed = entry.sequence[trim[0]]

    return trimmed


def main(args):
    """For each sequence, trim and write to output."""
    trimmed = []

    for entry in pysam.FastxFile(args.sequence):

        if 'trim=' not in entry.comment:
            continue

        trimmed.append([entry.name, trim(entry)])

    if not trimmed:
        return

    handler = get_output_handler(args.output)
    for name, seq in trimmed:
        handler.write(f">{name}\n{seq}\n")

    handler.close()


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("trim")
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
    args = argparse().parse_args()
    main(args)
