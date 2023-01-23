# Test data and scripts

This directory contains scripts and data to test the pipeline.
Some of the input files are large, so not all of them are tracked by git.
This means cloning the repository from Github will not download all the needed files.
Available testing scripts are:

## run_small.sh

Runs 2 samples with 10000 reads each.
These files are committed to the repositiory, so this command should always work
To use: run "test/run_small.sh" from the root of the repository.

## run_medium.sh

Runs 2 samples with all reads.
The reads are not committed to the repositiory.
To use: run "test/run_medium.sh" from the root of the repository.

