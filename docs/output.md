# nf-core/pathogensurveillance: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished.
All paths are relative to the top-level results directory which is located and named using the `--outdir` parameter.
Here, we will assume it is called `outdir`.


## Pipeline overview

The `pathogensurveillance` pipeline has many steps and each will produce output files for most of the steps.
Not all steps will be run for all input datasets.
For example, core gene phylogenies are only made for prokaryotes and busco gene phylogenies are only made for eukaryotes, so a dataset without both prokaryotes and eukaryotes will not have both of these outputs.
Below is the directory structure of the all possible outputs.

```
outdir
├── aligned_genes
│   ├── busco_genes
│   └── core_genes
├── aligned_reads
├── annotations
│   └── bakta
├── assemblies
│   ├── flye
│   └── spades
├── busco
├── downloads
│   ├── annotations
│   ├── assemblies
│   ├── databases
│   │   ├── bakta
│   │   └── busco
│   └── reads
├── fastp
├── metadata
├── pipeline_info
├── pirate
├── pocp
├── quality_control
│   ├── fastqc
│   ├── multiqc
│   ├── nanoplot
│   └── quast
├── reference_data
│   ├── considered
│   ├── downloaded
│   ├── indexes
│   │   ├── bgzip
│   │   ├── bwa
│   │   ├── faidx
│   │   ├── picard
│   │   └── tabix
│   └── selected
├── reports
├── report_group_data
├── sendsketch
├── sketch_comparisons
│   ├── ani_matricies
│   └── sketches
├── trees
│   ├── busco
│   ├── core
│   └── snp
└── variants
```

Within each of these directories there are files or directories named by the relevant ID type for the output.
For example, assemblies are named by the sample ID and reference indexes are named by the reference ID.
These IDs match those in metadata table in `metadata/sample_metadata.tsv` and `metadata/reference_metadata.tsv`, making it easy to automate downstream analysis of the data.
Additionally the [PathoSurveilR package](https://github.com/grunwaldlab/PathoSurveilR) can be used to automatically find and parse various output files given a top-level output directory (i.e. `outdir` in this example) for use in R.

Below is a more detailed description of each output directory.


### Aligned genes (`mafft`)

<details markdown="1">
<summary>Output files</summary>

- `aligned_genes/`
  - `busco_genes/`
    - `<gene ID>_aligned.fas`: FASTA files of aligned genes used in the BUSCO gene phylogenies.
  - `core_genes/`
    - `<gene ID>_aligned.fas`: FASTA files of aligned genes used in the core gene phylogenies.

</details>

FASTA files for each gene extracted from assemblies and aligned.
Contains sequences for both samples and references.


### Aligned reads (`bwa mem`)

<details markdown="1">
<summary>Output files</summary>

- `aligned_reads/`
  - `<Reference ID>_<Sample ID>.bam`: Alignments of reads to references in the BAM format.
  - `<Reference ID>_<Sample ID>.formatted.bam`: Quality filtered BAM files produced by `picard`.
  - `<Reference ID>_<Sample ID>.formatted.bam.csi`: Index for the above file produced by `samtools index`
  - `<Reference ID>_<Sample ID>.formatted.MarkDuplicates.metrics.txt`: Output from `picard MarkDuplicates`


</details>

Reads are aligned to references as part of the variant calling process used to compare samples with high resolution.
These read alignments are then filtered for quality and reformatted before being used to call variants.


### Prokaryotic gene annotations (Bakta)

<details markdown="1">
<summary>Output files</summary>

- `annotations/bakta/`
  - `<samplename>.gff3`: Annotations and sequences in GFF3 format
  - `<samplename>.gbff`: Annotations and sequences in (multi) GenBank format
  - `<samplename>.ffn`: Feature nucleotide sequences as FASTA
  - `<samplename>.fna`: Replicon/contig DNA sequences as FASTA
  - `<samplename>.embl`: Annotations and sequences in (multi) EMBL format
  - `<samplename>.faa`: CDS/sORF amino acid sequences as FASTA
  - `<samplename>_hypothetical.faa`: Further information on hypothetical protein CDS as simple human readable tab separated values
  - `<samplename>_hypothetical.tsv`: Hypothetical protein CDS amino acid sequences as FASTA
  - `<samplename>.tsv`: Annotations as simple human readble TSV
  - `<samplename>.txt`: Summary in TXT format

</details>

[Bakta](https://github.com/oschwengers/bakta) is a tool for the rapid and standardised annotation of bacterial genomes and plasmids from both isolates and MAGs.
It is used to annotate prokaryotic genomes for use in the core gene phylogeny.


### Assemblies (Spades and Flye)

<details markdown="1">
<summary>Output files</summary>

- `assemblies/`
  - `spades/`
    - `<Sample ID>.scaffolds.fa.gz`: Compressed assembled scaffolds in fasta format
    - `<Sample ID>.assembly.gfa.gz`: Compressed assembly graph in gfa format
    - `<Sample ID>.contigs.fa.gz`: Compressed assembled contigs in fasta forma
    - `<Sample ID>.spades.log`: Log file produced by `spades`
    - `<Sample ID>_filtered.fasta`: Quality filtered spades assembly
  - `flye/`
    - `<Sample ID>.assembly.fasta.gz`: Assembly in gzipped fasta format
    - `<Sample ID>.assembly_graph.gfa.gz`: Assembly graph in gzipped gfa format
    - `<Sample ID>.assembly_graph.gv.gz`: Assembly graph in gzipped gv format
    - `<Sample ID>.assembly_info.txt`: Information on the assembly
    - `<Sample ID>.flye.log`: Flye log file
    - `<Sample ID>.params.json`: Parameters used when running flye

</details>

These directories contain the output of whole genome assembly of samples using `spades` for short reads and `flye` for long reads.


### BUSCO

<details markdown="1">
<summary>output files</summary>

- `busco/`
  - `short_summary.specific.<busco_db>.<species_name>.fasta.txt`: completeness report in tsv format
  - `<species_name>-<busco_db>-busco.batch_summary.txt`: summarized completeness report in tsv format
  - `<sample id>-<database lineage>-busco`: directory with other busco results

</details>

BUSCO is used to extract genes from eukaryotic assemblies for phylogenetic analysis and assess assembly completeness.


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

</details>

This directory contains anything the pipeline downloads, such as assemblies, reads, and databases.


### Adapter trimming and quality control (`fastp`)

<details markdown="1">
<summary>output files</summary>

- `fastp/`
  - `<Sample ID>.fastp.fastq.gz`: Adapter trimmed FASTQ files
  - `<Sample ID>.fastp.html`: FASTP report
  - `<Sample ID>.fastp.json`: JSON data for the above report
  - `<Sample ID>.fastp.log`: Runtime log for FASTP

</details>

`fastp` is used to trim adapters and for other quality control.
It also produces a useful report on the quality of the sample.


### Sample and reference metadata

<details markdown="1">
<summary>output files</summary>

- `metadata/`
  - `sample_metadata.tsv`: A table with cleaned user sample metadata
  - `ref_metadata.tsv`: A table with cleaned user sample metadata

</details>

These files are the parsed and cleaned versions of the input data.
The IDs present in these tables are those used throughout the pipeline and might be different than what the user provided if they needed to be renamed to be compatible with use in file names.
These versions of the metadata should be used to automate any downstream analysis rather than the input metadata provided by the user.


### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides various reports relevant to the running and execution of the pipeline.
This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.


### Pirate

<details markdown="1">
<summary>output files</summary>

- `pirate/`
  - `<Report ID>_results`: Pirate output

</details>

Pirate is used to identify orthologous gene clusters, which are used later in the pipeline to create phylogenies of prokaryotes with the maximum number of shared genes without relying on annotations.


### Percentage of conserved proteins (POCP)

<details markdown="1">
<summary>Output files</summary>

- `pocp`
  - `<Report ID>_pocp.tsv`: A pairwise matrix of the POCP between all samples and references

</details>

POCP is calculated as a metric to compare samples to each other and to references in regards to shared gene content.


### Quality control reports

<details markdown="1">
<summary>output files</summary>

- `quality_control/`
  - `multiqc/`
    - `<Report ID>_multiqc>`: MultiQC outputs for samples in each report
  - `nanoplot/`: Nanoplot output reports and plots
  - `quast/`
    - `<Sample ID>`: Quast reports and associated data
  - `fastqc/`: FASTQC output reports

</details>

Various tools are used to check reads and assemblies for quality.
The outputs of these tools are compiled using MultiQC.


### Reference data

<details markdown="1">
<summary>Output files</summary>

- `reference_data/`
  - `considered/`
    - `<family>.json`: The metadata downloaded from the NCBI assembly database for all references considered
    - `<family>.tsv`: Select information from the above JSON file converted to a table for easy parsing
  - `downloaded/`
    - `<sample_id>.tsv`: The metadata for references selected for download
  - `selected/`
    - `<report group>_mapping_references.tsv`: The IDs of references used to align reads to during variant calling
    - `<report group>_core_references.tsv`: The IDs of references used to provide context in core gene phylogenies
    - `<report group>_busco_references.tsv`:  The IDs of references used to provide context in BUSCO gene phylogenies
  - `indexes/`
    - `bwa/`
      - `<Refernce ID>_bwa`: Index files used to align reads to references with `bwa mem`
    - `tabix/`
      - `<Report ID>_<Reference ID>.vcf.gz.tbi`: Index files created by `tabix`, which is part of samtools
    - `bgzip/`
      - `<Reference ID>.fasta.gz.gzi`: Index files created by `bgzip`, which is part of samtools
    - `faidx`:
      - `<Reference ID>.fasta.gz.fai`: Index files created by `faidx`, which is part of samtools
      - `<Reference ID>.fasta.gz.gzi`: Index files created by `faidx`, which is part of samtools
    - `picard`:
      - `<Reference ID>.fasta.dict`: Index files created by `picard CreateSequenceDictionary`.

</details>

The `pathogensurveillance` pipeline will select and download references automatically for use in multiple steps throughout the pipeline.
The `reference_data` folder contains information regarding references, including the metadata of all that were considered, the metadata of those downloaded, and the IDs of those selected for use in the analysis.


### Main reports

<details markdown="1">
<summary>output files</summary>

- `reports/`
  - `<Report ID>_report.html`: The primary output report of the pipeline

</details>

This is the primary output of the pipeline, containing the report meant to be understandable by non-bioinformaticians.


### Grouped report data

<details markdown="1">
<summary>output files</summary>

- `report_group_data/`
  - `<Report ID>_inputs`: A folder containing formatted outputs from the pipeline used in the main report.

</details>

This is the directory used to create the main report for each report group.
It contains selected and renamed outputs from the pipeline present in other output folders, but organized by report group.


### BBMap Sendsketch results

<details markdown="1">
<summary>Output files</summary>

- `sendsketch/`
  - `<Sample ID>.txt`: Table returned by BBmap `sendsketch` with initial identifications of samples.

</details>

Tables with information used to make initial identifications of samples from the BBMap `sendsketch` tool.


### Hash-based comparisons

<details markdown="1">
<summary>output files</summary>

- `sketch_comparisons/`
  - `ani_matricies/`
    - `<Report ID>_comp.csv`: ANI similarity matrix in CSV format made by `sourmash compare`
    - `<Report ID>_comp.npy`: ANI similarity matrix in NumPy format made by `sourmash compare`
    - `<Report ID>_comp.npy.labels.txt`: Labels for the above file made by `sourmash compare`
  - `sketches/`
    - `<Sample ID or Reference ID>.sig`: FracMinHash signature of the given sequence made by `sourmash sketch`


</details>

In order to select references to use with samples and provide a rough identification, all samples and references are sketched with `sourmash sketch` and all pairwise comparisons of sketches are made with `sourmash compare`.


### Trees (`iqtree2`)

<details markdown="1">
<summary>output files</summary>

- `trees/`
  - `busco/`
    - `<Report ID>_<Cluster ID>.treefile`: Tree in Newick format inferred from BUSCO genes by `iqtree2`
  - `core/`
    - `<Report ID>_<Cluster ID>.treefile`: Tree in Newick format inferred from core genes by `iqtree2`
  - `snp/`
    - `<Report ID>_<Cluster ID>.treefile`: Tree in Newick format inferred from variants by `iqtree2`

</details>

Various trees are produced by the pipeline to compare the samples to references and to each other.
To put samples in context of reference genomes and provide data that can be useful in identification, core genes from prokaryotes and BUSCO genes from eukaryotes are used to produce trees with `iqtree2`.
SNPs identified by variant calling are also used to create a tree with `iqtree2` for high-resolution sample comparison.


### Variants (`graphtyper genotype`)

<details markdown="1">
<summary>output files</summary>

- `variants/`
  - `<Report ID>_<Reference ID>.vcf.gz`: The variants for all samples aligned to this reference produced by `graphtyper genotype`
  - `<Report ID>_<Reference ID>.vcf.gz.tbi`: The index files for the variants
  - `<Report ID>_<Reference ID>variantfiltration.vcf.gz`: The filtered variants for all samples aligned to this reference produced by `graphtyper genotype`
  - `<Report ID>_<Reference ID>variantfiltration.vcf.gz.tbi`: The index files for the filtered variants
  - `<Report ID>_<Reference ID>.vcffilter.vcf.gz`:
  - `<Report ID>_<Reference ID>.fasta`: FASTA file with values for each variable site concatenated

</details>

Variants are called against selected references to do a high-resolution comparison of samples.


