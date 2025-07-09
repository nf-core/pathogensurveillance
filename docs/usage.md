# nf-core/pathogensurveillance: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/pathogensurveillance/usage](https://nf-co.re/pathogensurveillance/usage)

## Introduction

This pipeline is designed to be as simple to run as possible while still being highly customizable.
If you are running the pipeline for the first time, consider running one of the test profiles described in the Introduction page (i.e. the README on Github).

## Samplesheet input

The primary input to the pipeline is a TSV (tab-separated value) or CSV (comma-separated value) file, specified using the `--input` option.
This can be made in a spreadsheet program like LibreOffice Calc or Microsoft Excel by exporting to TSV/CSV.
Use this parameter to specify its location:

```bash
--input <PATH TO SAMPLESHEET FILE>
```

Columns can be in any order and unneeded columns can be left out or left blank.
Column names are case insensitive and spaces are equivalent to underscores and can be left out.
Only a single column containing either paths to raw sequence data, SRA (Sequence Read Archive) accessions, or NCBI queries to search the SRA is required and each sample can have values in different columns.
For example, the following is a valid input:

```csv title="samplesheet.csv"
sample_id,path_1,path_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,
```

Any columns not recognized by `pathogensurveillance` will be ignored, allowing users to adapt existing sample metadata table by adding new columns.
Below is a description of each column used by `pathogensurveillance`:

| Column             | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| ------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample_id`        | The unique identifier for each sample. This will be used in file names to distinguish samples in the output. Each sample ID must correspond to a single source of sequence data (e.g. the `path` and `ncbi_accession` columns), although the same sequence data can be used by different IDs. Any values supplied that correspond to different sources of sequence data or contain characters that cannot appear in file names will be modified automatically. If not supplied, it will be inferred from the `path`, `ncbi_accession`, or `name` columns. |
| `name`             | A human-readable label for the sample that is used in plots and tables. If not supplied, it will be inferred from `sample_id`.                                                                                                                                                                                                                                                                                                                                                                                                                            |
| `description`      | A longer human-readable label that is used in plots and tables. If not supplied, it will be inferred from `name`.                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| `path`             | Path to input sequence data, typically gzipped FASTQ files. When paired end sequencing is used, this is used for the forward read's data and `path_2` is used for the reverse reads. This can be a local file path or a URL to an online location. The `sequence_type` column must have a value.                                                                                                                                                                                                                                                          |
| `path_2`           | Path to the FASTQ files for the reverse read when paired-end sequencing is used. This can be a local file path or a URL to an online location. The `sequence_type` column must have a value.                                                                                                                                                                                                                                                                                                                                                              |
| `ncbi_accession`   | An SRA accession ID for reads to be downloaded and used as samples. Values in the `sequence_type` column will be looked up if not supplied.                                                                                                                                                                                                                                                                                                                                                                                                               |
| `ncbi_query`       | A valid NCBI search query to search the SRA for reads to download and use as samples. This will result in an unknown number of samples being analyzed. The total number downloaded is limited by the `ncbi_query_max` column. Values in the `sample_id`, `name`, and `description` columns will be append to that supplied by the user. Values in the `sequence_type` column will be looked up and does not need to be supplied by the user.                                                                                                              |
| `ncbi_query_max`   | The maximum number or percentage of samples downloaded for the corresponding query in the `ncbi_query` column. Adding a `%` to the end of a number indicates a percentage of the total number of results instead of a count. A random of subset of results will be downloaded if `ncbi_query_max` is less than "100%" or the total number of results.                                                                                                                                                                                                     |
| `sequence_type`    | The type of sequencing used to produce reads for the `reads_1` and `reads_2` columns. Valid values include anything containing the words "illumina", "nanopore", or "pacbio". Will be looked up automatically for `ncbi_accession` and `ncbi_query` inputs but must be supplied by the user for `path` inputs.                                                                                                                                                                                                                                            |
| `report_group_ids` | How to group samples into reports. For every unique value in this column a report will be generated. Samples can be assigned to multiple reports by separating group IDs by ";". For example `all;subset` will put the sample in both `all` and `subset` report groups. Samples will be added to a default group if this is not supplied.                                                                                                                                                                                                                 |
| `color_by`         | The names of other columns that contain values used to color samples in plots and figures in the report. Multiple column names can be separated by ";". Specified columns can contain either categorical factors or specific colors, specified as a hex code. By default, samples will be one color and references another.                                                                                                                                                                                                                               |
| `ploidy`           | The ploidy of the sample. Should be a number. Defaults to "1".                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| `enabled`          | Either "TRUE" or "FALSE", indicating whether the sample should be included in the analysis or not. Defaults to "TRUE".                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| `ref_group_ids`    | One or more reference group IDs separated by ";". These are used to supply specific references to specific samples. These IDs correspond to IDs listed in the `ref_group_ids` or `ref_id` columns of the reference metadata TSV.                                                                                                                                                                                                                                                                                                                          |

Additionally, users can supply a reference metadata TSV/CSV that can be used to assign custom references to particular samples using the `--reference_data` option.
If not provided, the pipeline will download and choose references to use automatically.
References are assigned to samples if they share a reference group ID in the `ref_group_ids` columns that can appear in both input TSVs/CSVs.
For example, the following file content given to `--reference_data` would assign all samples to a given reference based on an NCBI accession:

```csv title="samplesheet.csv"
ref_group_ids,ncbi_accession,ref_contextual_usage
my_refs,GCA_000513215.1,required
my_refs,GCF_039778255.1,required
```

To make this have an effect, the `--input` spreadsheet must have the `ref_group_ids` column defined with the same ID:

```csv title="samplesheet_with_ref.csv"
sample_id,path_1,path_2,ref_group_ids
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,my_refs
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz,my_refs
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz,my_refs
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,,my_refs
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,,my_refs
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,,my_refs
```

The reference metadata TSV or the sample metadata TSV can have the following columns:

| Column                 | Description                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| ---------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `ref_group_ids`        | One or more reference group IDs separated by ";". These are used to group references and supply an ID that can be used in the `ref_group_ids` column of the sample metadata TSV/CSV to assign references to particular samples.                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| `ref_id`               | The unique identifier for each user-defined reference genome. This will be used in file names to distinguish samples in the output. Each reference ID must correspond to a single source of reference data (The `ref_path`, `ref_ncbi_accession`, and `ref_ncbi_query` columns), although the same reference data can be used by multiple IDs. Any values that correspond to different sources of reference data or contain characters that cannot appear in file names will be modified automatically. If not supplied, it will be inferred from the `path`, `ref_name` columns or supplied automatically when `ref_ncbi_accession` or `ref_ncbi_query` are used. |
| `ref_name`             | A human-readable label for user-defined reference genomes that is used in plots and tables. If not supplied, it will be inferred from `ref_id`. It will be supplied automatically when the `ref_ncbi_query` column is used.                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| `ref_description`      | A longer human-readable label for user-defined reference genomes that is used in plots and tables. If not supplied, it will be inferred from `ref_name`. It will be supplied automatically when the `ref_ncbi_query` column is used.                                                                                                                                                                                                                                                                                                                                                                                                                               |
| `ref_path`             | Path to user-defined reference genomes for each sample. This can be a local file path or a URL to an online location.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| `ref_ncbi_accession`   | RefSeq accession ID for a user-defined reference genome. These will be automatically downloaded and used as input.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| `ref_ncbi_query`       | A valid NCBI search query to search the assembly database for genomes to download and use as references. This will result in an unknown number of references being downloaded. The total number downloaded is limited by the `ref_ncbi_query_max` column. Values in the `ref_id`, `ref_name`, and `ref_description` columns will be append to that supplied by the user.                                                                                                                                                                                                                                                                                           |
| `ref_ncbi_query_max`   | The maximum number or percentage of references downloaded for the corresponding query in the `ref_ncbi_query` column. Adding a `%` to the end of a number indicates a percentage of the total number of results instead of a count. A random of subset of results will be downloaded if `ncbi_query_max` is less than "100%" or the total number of results.                                                                                                                                                                                                                                                                                                       |
| `ref_primary_usage`    | Controls how the reference is used in the analysis in cases where a single "best" reference is required, such as for variant calling. Can be one of "optional" (can be used if selected by the analysis), "required" (will always be used), "exclusive" (only those marked "exclusive" will be used), or "excluded" (will not be used).                                                                                                                                                                                                                                                                                                                            |
| `ref_contextual_usage` | Controls how the reference is used in the analysis in cases where multiple references are required to provide context for the samples, such as for phylogeny. Can be one of "optional" (can be used if selected by the analysis), "required" (will always be used), "exclusive" (only those marked "exclusive" will be used), or "excluded" (will not be used).                                                                                                                                                                                                                                                                                                    |
| `ref_color_by`         | The names of other columns that contain values used to color references in plots and figures in the report. Multiple column names can be separated by ";". Specified columns can contain either categorical factors or specific colors, specified as a hex code. By default, samples will be one color and references another.                                                                                                                                                                                                                                                                                                                                     |
| `ref_enabled`          | Either `TRUE` or `FALSE`, indicating whether the reference should be included in the analysis or not. Defaults to `TRUE`.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |

## Running the pipeline

> - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile docker`.
> - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
> - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
> - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/pathogensurveillance -profile <REPLACE WITH PACKAGE MANAGER> -resume --input <REPLACE WITH TSV/CSV> --outdir <REPLACE WITH OUTPUT PATH>
```

Where:

- `<REPLACE WITH PACKAGE MANAGER>` is one of docker, singularity, podman, shifter, charliecloud, or conda
- `<REPLACE WITH TSV/CSV>` is the path to the input samplesheet
- `<REPLACE WITH OUTPUT PATH>` is the path to where to save the output

An actual command might look like this:

```bash
nextflow run nf-core/pathogensurveillance -profile docker -resume --input ./sample_metadata.tsv --outdir ./results
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
path_surveil_data   # Where download reads and references are stored for reuse. Can be changed with the `--data_dir` parameter
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.
Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/pathogensurveillance -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version.
When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since
To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/pathogensurveillance
```

### API keys

The pipeline uses multiple tools that communicate with external services, such as NCBI.
These services impose usage limits in order to ensure fair access to services.
The pipeline is designed to not exceed these rate limits.
Users can register for free API keys to increase these rate limits and the pipeline with automatically adjust to reflect this.
Currently, only the NCBI API key is used by the pipeline.
It can be set using [Nextflow secrets](https://www.nextflow.io/docs/latest/secrets.html):

```bash
nextflow secrets set NCBI_API_KEY INSERT_YOUR_KEY_HERE
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/pathogensurveillance releases page](https://github.com/nf-core/pathogensurveillance/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
