# Test data and scripts

This directory contains data to test the pipeline.
Some of the input files are large, so not all of them are tracked by git.
This means cloning the repository from Github will not download all the needed files.
Available testing scripts are:

## small data set

Reads used are in "data/reads_small"

3 samples with 10000 reads each.
These files are committed to the repositiory, so this command should always work
Use the following command from the root of the repository to run this test dataset:

nextflow run main.nf -profile test,docker

The output will be in "test/output"
