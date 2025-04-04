# nf-core/pathogensurveillance: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The `pathogensurveillance` pipeline has many steps and each will produce output files in a directory named by the step.
Not all steps will be run for all input datasets.
Below is a list of the most important outputs created.


### Prokaryotic gene annotations (Bakta)

<details markdown="1">
<summary>Output files</summary>

- `annotations/bakta/`
  - `<samplename>.gff3`: annotations & sequences in GFF3 format
  - `<samplename>.gbff`: annotations & sequences in (multi) GenBank format
  - `<samplename>.ffn`: feature nucleotide sequences as FASTA
  - `<samplename>.fna`: replicon/contig DNA sequences as FASTA
  - `<samplename>.embl`: annotations & sequences in (multi) EMBL format
  - `<samplename>.faa`: CDS/sORF amino acid sequences as FASTA
  - `<samplename>_hypothetical.faa`: further information on hypothetical protein CDS as simple human readble tab separated values
  - `<samplename>_hypothetical.tsv`: hypothetical protein CDS amino acid sequences as FASTA
  - `<samplename>.tsv`: annotations as simple human readble TSV
  - `<samplename>.txt`: summary in TXT format

> Descriptions taken from the [Bakta documentation](https://github.com/oschwengers/bakta#output).

</details>

[Bakta](https://github.com/oschwengers/bakta) is a tool for the rapid & standardised annotation of bacterial genomes and plasmids from both isolates and MAGs. It provides dbxref-rich, sORF-including and taxon-independent annotations in machine-readable JSON & bioinformatics standard file formats for automated downstream analysis.


### Reference data

<details markdown="1">
<summary>Output files</summary>

- `reference_data/`
  - `considered/`
    - `<family>.json`: The metadata downloaded from the NCBI assembly database for all references considered
    - `<family>.tsv`: Select information from the above JSON file converted to a table
  - `downloaded/`
    - `<sample_id>.tsv`: The metadata for references selected for download
  - `selected/`
    - `<report group>_mapping_references.tsv`: The IDs of references used to align reads to during variant calling
    - `<report group>_core_references.tsv`: The IDs of references used to provide context in core gene phylogenies
    - `<report group>_busco_references.tsv`:  The IDs of references used to provide context in BUSCO gene phylogenies

> Descriptions taken from the [Bakta documentation](https://github.com/oschwengers/bakta#output).

</details>

The `pathogensurveillance` pipeline will select and download references automatically for use in multiple steps throughout the pipeline.
The `reference_data` folder contains information regarding references, including the metadata of all that were considered, the metada of those downloaded, and the IDs of those selected for use in the analysis.


### Downloads

<details markdown="1">
<summary>Output files</summary>

- `downloads/`
  - `assemblies/`
    - `<reference ID>.fasta.gz`: FASTA files of assemblies
  - `annotations/`
    - `<reference ID>.gff.gz`: GFF files of annotations corresponding to assemblies
  - `reads/`
    - `<sample ID>.gff.gz`: FASTQ files of reads
  - `databases/`
    - `bakta/`
      - `db-<size>`: The database used for annotation of assemblies using Bakta. The `<size>` can be `full` or `light`
    - `busco/`
      - `busco_downloads`: The database used for identification of single copy orthologs using BUSCO

> Descriptions taken from the [Bakta documentation](https://github.com/oschwengers/bakta#output).

</details>

This directory contains anything the pipeline downloads, such as assemblies, reads, and databases.









### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
