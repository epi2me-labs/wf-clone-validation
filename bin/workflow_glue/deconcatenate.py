#!/usr/bin/env python
"""Go deconcatenate your sequences."""
import sys

import mappy as mp
from .util import wf_parser  # noqa: ABS101


def get_output_handler(path):
    """Open file path or stdout."""
    if path != '-':
        fh = open(path, 'w')
    else:
        fh = sys.stdout
    return fh


def get_aligner(reference):
    """Return a ready aligner."""
    aligner = mp.Aligner(seq=reference, preset='asm5')
    if not aligner:
        raise Exception("ERROR: failed to load/build index")
    return aligner


def align_self(seq):
    """Split read and align one half to the other."""
    half = len(seq) // 2
    first, second = seq[0:half], seq[half:]

    aligner = get_aligner(first)
    hits = [hit for hit in aligner.map(second)]

    return hits, first, second


def deconcatenate(seq):
    """Self-align to remove duplicate regions."""
    finished = False
    iteration = 0
    trimmed_assm = seq

    while not finished:
        iteration += 1
        sys.stdout.write(f"Trimming sequence... Round {iteration}")
        hits, first, second = align_self(trimmed_assm)

        if len(hits) == 1:
            sys.stdout.write("> Single self-alignment detected.")

        elif len(hits) > 1:
            sys.stdout.write("> Multiple self-alignments detected.")
            # Tested variations of this, but if works...
            hits = [hit for hit in hits if hit.q_st < 5]
            if not hits:
                sys.stdout.write(
                    "> No self-alignments match criteria, stopping here.")
                finished = True

        else:
            sys.stdout.write("> No self-alignments, stopping here.")
            finished = True
            break
        try:
            hit = hits[0]
            trimmed_assm = second[:hit.q_en] + first[hit.r_en:]
        except IndexError:
            trimmed_assm = seq

    return trimmed_assm


def main(args):
    """For each sequence, deconcatenate and write to output."""
    sequence_fasta = args.sequence
    output = args.output
    corrected = []
    for name, seq, _ in mp.fastx_read(sequence_fasta):
        corrected.append([name, deconcatenate(seq)])

    if not corrected:
        return

    handler = get_output_handler(output)
    for n, s in corrected:
        handler.write(f">{n}\n{s}\n")
    handler.close()


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("deconcatenate")
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


if __name__ == '__main__':
    args = argparser().parse_args()
    main(args)
