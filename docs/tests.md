## Overview

This pipeline supports several full-pipeline tests through nf-test. All run a complete pipeline, but many are intended to test specific features during development.

## Setup

- Follow nf-test documentation to install the nf-test executable: https://www.nf-test.com/installation/.
- If needed, edit nf-test.config to configure nf-test-specific settings (testing configuration file for Nextflow, default profile used, plugins, etc.).
- Also edit tests/nextflow.config (default Nextflow configuration file for testing) to configure extra Nextflow profiles and insititutional executor settings.

## Usage

### Run Tests

Run all tests as such:

`nf-test test`

Or, to run a specific test,

`nf-test test tests/my_test`

### Cache Clearing and Snapshot Updates

To clear nf-test cache (nf-test-work directory) use:

`nf-test clean`

To run tests and regenerate their snapshots, run:

`nf-test tests/my_test --update-snapshot`

### Other Arguments

- `--profile +my_profile`: Specify additional Nextflow profiles (`+my_profile` will be additive to the default profile, `my_profile` will replace the default profile).
- `--verbose`: View active pipeline execution steps.
- `--debug`: Print debug messages.
- `--tag`: Subset test runs by category

## Pipeline Tags

Most tests for this pipeline have tags related to the pipeline section or feature they test. Apply one or more to the test command to modify the set of tests being run.

- `input`: Tests an input component of the pipeline, skips core phylogeny.
- `subworkflow`: Tests a specific pipeline subworkflow.
- `genomic`: Tests a specific genomic category of samples.
- `full_pipeline`: Guaranteed to run complete pipeline logic.
- `dev`: Small subset of tests intended for minimal but comprehensive testing during pipeline development.
