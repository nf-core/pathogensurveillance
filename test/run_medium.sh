#! /usr/bin/env bash

# This command runs the pipeline using the test data included in the test/data directory
# NOTE: this assumes that the command is run from the root of the pipeline repositiory,
#     not in the test directory
nextflow run main.nf -profile docker --input test/data/reads_medium/metadata.csv --fasta test/data/genome.fa --outdir test/output
