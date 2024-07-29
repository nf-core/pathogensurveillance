<!--
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-pathogensurveillance_logo_dark.png">
  <source media="(prefers-color-scheme: light)" srcset="docs/images/nf-core-pathogensurveillance_logo_light.png">
  <img alt="nf-core/pathogensurveillance" src="docs/images/nf-core-pathogensurveillance_logo_light.png">
</picture>
-->

[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/pathogensurveillance/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/pathogensurveillance)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23pathogensurveillance-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/pathogensurveillance)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## NOTE: THIS PROJECT IS UNDER DEVELOPMENT AND MAY NOT FUNCTION AS EXPECTED UNTIL THIS MESSAGE GOES AWAY

## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->

**nf-core/pathogensurveillance** is a population genomic pipeline for pathogen diagnosis, variant detection, and biosurveillance.
The pipeline accepts the paths to raw reads for one or more organisms (in the form of a CSV file) and creates reports in the form of interactive HTML reports or PDF documents.
Significant features include the ability to analyze unidentified eukaryotic and prokaryotic samples, creation of reports for multiple user-defined groupings of samples, automated discovery and downloading of reference assemblies from NCBI RefSeq, and rapid initial identification based on k-mer sketches followed by a more robust core genome phylogeny and SNP-based phylogeny.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner.
It uses Docker/Singularity containers making installation trivial and results highly reproducible.
The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies.
Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure.
This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world data sets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources.The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/pathogensurveillance/results).

## Pipeline summary

![Pipeline flowchart](docs/posters/pipeline_diagram.png)

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```bash
   nextflow run nf-core/pathogensurveillance -r dev -profile RUN_TOOL,xanthomonas_small -resume --out_dir test_output
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`RUN_TOOL` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile xanthomonas_small,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

   ```bash
   nextflow run nf-core/pathogensurveillance -r dev -profile RUN_TOOL -resume --sample_data <CSV> --out_dir <OUTDIR> --download_bakta_db
   ```

```



## Documentation

The nf-core/pathogensurveillance pipeline comes with documentation about the pipeline [usage](https://nf-co.re/pathogensurveillance/usage), [parameters](https://nf-co.re/pathogensurveillance/parameters) and [output](https://nf-co.re/pathogensurveillance/output).

### Input format

The primary input to the pipeline is a CSV (comma comma-separated value) file, specified using the `--sample_data` option.
This can be made in a spreadsheet program like LibreOffice Calc or Microsoft Excel by exporting to CSV.
Columns can be in any order and unneeded columns can be left out or left blank.
Column names are case insensitive and spaces are equivalent to underscores.
Only a single column containing either paths to raw sequence data, SRA (Sequence Read Archive) accessions, or NCBI queries to search the SRA is required and each sample can have values in different columns.
Any columns not recognized by `pathogensurveillance` will be ignored, allowing users to adapt existing sample metadata table by adding new columns.
Below is a description of each column used by `pathogensurveillance`:

* **sample_id**: The unique identifier for each sample. This will be used in file names to distinguish samples in the output. Each sample ID must correspond to a single source of sequence data (e.g. the `path` and `ncbi_accession` columns), although the same sequence data can be used by different IDs. Any values supplied that correspond to different sources of sequence data or contain characters that cannot appear in file names (\/:*?"<>| .) will be modified automatically. If not supplied, it will be inferred from the `path`, `ncbi_accession`, or `name` columns.
* **name**: A human-readable label for the sample that is used in plots and tables. If not supplied, it will be inferred from `sample_id`.
* **path**: Path to input sequence data, typically gzipped FASTQ files. When paired end sequencing is used, this is used for the forward read's data and `path_2` is used for the reverse reads. This can be a local file path or a URL to an online location. The `sequence_type` column must have a value.
* **path_2**: Path to the FASTQ files for the reverse read when paired-end sequencing is used. This can be a local file path or a URL to an online location. The `sequence_type` column must have a value.
* **ncbi_accession**: An SRA accession ID for reads to be downloaded and used as samples. Values in the `sequence_type` column will be looked up if not supplied.
* **ncbi_query**: A valid NCBI search query to search the SRA for reads to download and use as samples. This will result in an unknown number of samples being analyzed. The total number downloaded is limited by the `ncbi_query_max` column.  Values in the `sample_id`, `name`, and `description` columns will be append to that supplied by the user. Values in the `sequence_type` column will be looked up and does not need to be supplied by the user.
* **ncbi_query_max**: The maximum number or percentage of samples downloaded for the corresponding query in the `ncbi_query` column. Adding a `%` to the end of a number indicates a percentage of the total number of results instead of a count. A random of subset of results will be downloaded if `ncbi_query_max` is less than "100%" or the total number of results.
* **sequence_type**: The type of sequencing used to produce reads for the `reads_1` and `reads_2` columns. Valid values include anything containing the words "illumina", "nanopore", or "pacbio". Will be looked up automatically for `ncbi_accession` and `ncbi_query` inputs but must be supplied by the user for `path` inputs.
* **report_group_ids**: How to group samples into reports. For every unique value in this column a report will be generated. Samples can be assigned to multiple reports by separating group IDs by ";". For example `all;subset` will put the sample in both `all` and `subset` report groups. Samples will be added to a default group if this is not supplied.
* **color_by**: The names of other columns that contain values used to color samples in plots and figures in the report. Multiple column names can be separated by ";". Specified columns can contain either categorical factors or specific colors, specified as a hex code. By default, samples will be one color and references another.
* **ploidy**: The ploidy of the sample. Should be a number. Defaults to "1".
* **enabled**: Either "TRUE" or "FALSE", indicating whether the sample should be included in the analysis or not. Defaults to "TRUE".
* **ref_group_ids**: One or more reference group IDs separated by ";". These are used to supply specific references to specific samples. These IDs correspond to IDs listed in the `ref_group_ids` or `ref_id` columns of the reference metadata CSV.

Additionally, users can supply a reference metadata CSV that can be used to assign custom references to particular samples using the `--reference_data` option.
If not provided, the pipeline will download and choose references to use automatically.
References are assigned to samples if they share a reference group ID in the `ref_group_ids` columns that can appear in both input CSVs.
The reference metadata CSV or the sample metadata CSV can have the following columns:

* **ref_group_ids**: One or more reference group IDs separated by ";". These are used to group references and supply an ID that can be used in the `ref_group_ids` column of the sample metadata CSV to assign references to particular samples. * **ref_id**: The unique identifier for each user-defined reference genome. This will be used in file names to distinguish samples in the output. Each reference ID must correspond to a single source of reference data (The `ref_path`, `ref_ncbi_accession`, and `ref_ncbi_query` columns), although the same reference data can be used by multiple IDs. Any values that correspond to different sources of reference data or contain characters that cannot appear in file names (\/:*?"<>| .) will be modified automatically. If not supplied, it will be inferred from the `path`, `ref_name` columns or supplied automatically when `ref_ncbi_accession` or `ref_ncbi_query` are used.
* ref_id: The unique identify for each reference input. This will be used in file names to distinguish references in the output. Each sample ID must correspond to a single source of reference data (e.g. the `ref_path` and `ref_ncbi_accession` columns), although the same sequence data can be used by different IDs. Any values supplied that correspond to different sources of reference data or contain characters that cannot appear in file names (\/:*?"<>| .) will be modified automatically. If not supplied, it will be inferred from the `ref_path`, `ref_ncbi_accession`, or `ref_name` columns.
* **ref_name**: A human-readable label for user-defined reference genomes that is used in plots and tables. If not supplied, it will be inferred from `ref_id`. It will be supplied automatically when the `ref_ncbi_query` column is used.
* **ref_description**: A longer human-readable label for user-defined reference genomes that is used in plots and tables. If not supplied, it will be inferred from `ref_name`. It will be supplied automatically when the `ref_ncbi_query` column is used.
* **ref_path**: Path to user-defined reference genomes for each sample. This can be a local file path or a URL to an online location.
* **ref_ncbi_accession**: RefSeq accession ID for a user-defined reference genome. These will be automatically downloaded and used as input.
* **ref_ncbi_query**: A valid NCBI search query to search the assembly database for genomes to download and use as references. This will result in an unknown number of references being downloaded. The total number downloaded is limited by the `ref_ncbi_query_max` column.  Values in the `ref_id`, `ref_name`, and `ref_description` columns will be append to that supplied by the user.
* **ref_ncbi_query_max**: The maximum number or percentage of references downloaded for the corresponding query in the `ref_ncbi_query` column. Adding a `%` to the end of a number indicates a percentage of the total number of results instead of a count. A random of subset of results will be downloaded if `ncbi_query_max` is less than "100%" or the total number of results.
* **ref_primary_usage**: Controls how the reference is used in the analysis in cases where a single "best" reference is required, such as for variant calling. Can be one of "optional" (can be used if selected by the analysis), "required" (will always be used), "exclusive" (only those marked "exclusive" will be used), or "excluded" (will not be used).
* **ref_contextual_usage**: Controls how the reference is used in the analysis in cases where multiple references are required to provide context for the samples, such as for phylogeny. Can be one of "optional" (can be used if selected by the analysis), "required" (will always be used), "exclusive" (only those marked "exclusive" will be used), or "excluded" (will not be used).
* **ref_color_by**: The names of other columns that contain values used to color references in plots and figures in the report. Multiple column names can be separated by ";". Specified columns can contain either categorical factors or specific colors, specified as a hex code. By default, samples will be one color and references another.
* **ref_enabled**: Either "TRUE" or "FALSE", indicating whether the reference should be included in the analysis or not. Defaults to "TRUE".

## Credits

nf-core/pathogensurveillance was originally written by Zachary S.L. Foster, Camilo Parada-Rojas, Martha Sudermann, Nicholas C. Cauldron, Fernanda I. Bocardo, Ricardo Alcalá-Briseño, Hung Phan, Jeﬀ H. Chang, Niklaus J. Grünwald.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#pathogensurveillance` channel](https://nfcore.slack.com/channels/pathogensurveillance) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/pathogensurveillance for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
