#!/usr/bin/env python


"""Extract the relevant taxonomic information from the result of bbmap/sendsketch."""


import argparse
import csv
import sys
import logging
from pathlib import Path

logger = logging.getLogger()


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract taxonomic information from sendsketch output to determine which pipeline to use.",
        epilog="Example: python initial_classification.py sendsketch_result.txt",
    )
    parser.add_argument(
        "input",
        type=Path,
        help="Sendsketch output in TSV format.",
    )
    parser.add_argument(
        "taxon_out",
        type=Path,
        help="The file to write the taxon name to.",
    )
    parser.add_argument(
        "class_out",
        type=Path,
        help="The file to write the classification heirarchy to.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def extract_taxon_name(input_path, out_path="taxon.txt"):
    """Get coarse taxonomic classification from sendsketch output and save to file"""
    with open(input_path) as handle:
        for _ in range(2):  # Skip first 2 lines 
            next(handle)
        tsv_reader = csv.DictReader(handle, delimiter='\t')
        first_line = next(tsv_reader)
        taxon = first_line['taxonomy'].split(';')[0].split(':')[1]
    with open(out_path, "w") as handle:
        handle.write(taxon)

def extract_classification(input_path, out_path="class.txt"):
    """Get coarse taxonomic classification from sendsketch output and save to file"""
    with open(input_path) as handle:
        for _ in range(2):  # Skip first 2 lines 
            next(handle)
        tsv_reader = csv.DictReader(handle, delimiter='\t')
        first_line = next(tsv_reader)
        output = first_line['taxonomy']
    with open(out_path, "w") as handle:
        handle.write(output)

def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.input.is_file():
        logger.error(f"The given input file {args.input} was not found!")
        sys.exit(2)
    extract_taxon_name(args.input, args.taxon_out)
    extract_classification(args.input, args.class_out)


if __name__ == "__main__":
    sys.exit(main())

