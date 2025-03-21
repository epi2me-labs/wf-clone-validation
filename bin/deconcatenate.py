#!/usr/bin/env python
"""Go deconcatenate your sequences."""
import argparse
import sys

import mappy as mp


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


def deconcatenate(seq, approx_size):
    """Self-align to remove duplicate regions."""
    finished = False
    iteration = 0
    trimmed_assm = seq
    while not finished:
        iteration += 1
        sys.stdout.write(f"Trimming sequence... Round {iteration}\n")
        approx_size = int(approx_size)
        upper_limit = approx_size * 1.2
        lower_limit = approx_size * 0.8
        hits, first, second = align_self(trimmed_assm)
        if lower_limit < len(trimmed_assm) < upper_limit:
            sys.stdout.write(
                "Approx size is as expected, stopping here.\n")
            finished = True
            break
        elif len(hits) == 1:
            sys.stdout.write("> Single self-alignment detected.\n")
        elif len(hits) > 1:
            sys.stdout.write("> Multiple self-alignments detected.\n")
            # Tested variations of this, but if works...
            hits = [hit for hit in hits if hit.q_st < 5]
            if not hits:
                sys.stdout.write(
                    "> No self-alignments match criteria, stopping here.\n")
                finished = True
                break

        else:
            sys.stdout.write("> No self-alignments, stopping here.\n")
            finished = True
            break
        hit = hits[0]
        if hit.r_st < hit.q_st:
            trimmed_assm = second[:hit.q_en] + first[hit.r_en:]
        else:
            trimmed_assm = first[:hit.r_en] + second[hit.q_en:]

    return trimmed_assm


def main(args):
    """For each sequence, deconcatenate and write to output."""
    sequence_fasta = args.sequence
    output = args.output
    corrected = []
    for name, seq, _ in mp.fastx_read(sequence_fasta):
        corrected.append([name, deconcatenate(seq, args.approx_size)])

    if not corrected:
        return

    handler = get_output_handler(output)
    for n, s in corrected:
        handler.write(f">{n}\n{s}\n")
    handler.close()


def argparser():
    """Argument parser for entrypoint."""
    parser = argparse.ArgumentParser(
        "deconcatenate",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False)
    parser.add_argument(
        dest="sequence",
        help="File in .FASTA format containing a single sequence/contig."
    )
    parser.add_argument(
        "--approx_size",
        dest="approx_size",
        help="Approx plasmid size."
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
