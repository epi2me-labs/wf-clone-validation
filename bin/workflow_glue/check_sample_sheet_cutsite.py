#!/usr/bin/env python
"""Check if a sample sheet with cutsite is valid."""
import csv
import sys

from .util import get_named_logger, wf_parser  # noqa: ABS101


def main(args):
    """Run the entry point."""
    logger = get_named_logger("checkSheetCondition")
    with open(args.sample_sheet, "r") as f:
        csv_reader = csv.DictReader(f)
        for i in csv_reader:
            if any(x not in i['cut_site'] for x in 'AGCT'):
                sys.exit(
                    "Cut site column must not contain base pairs "
                    "other than AGCT.")
            if 5 > len(i['cut_site']) > 30:
                sys.exit(
                    f"Cut sites must all be between 5-30 bp long.{len(i['cut_site'])}")
    logger.info(f"Checked sample sheet for cutsite column {args.sample_sheet}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("check_sample_sheet_cutsite")
    parser.add_argument("sample_sheet", help="Sample sheet to check")
    return parser
