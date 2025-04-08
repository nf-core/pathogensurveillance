# nf-core/pathogensurveillance: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The `pathogensurveillance` pipeline has many steps and each will produce output files in a directory named by the step.
Not all steps will be run for all input datasets.
Below is a list of the most important outputs created.

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
├── fastqc
├── main_report_input
├── pipeline_info
├── pirate
├── pocp
├── quality_control
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
```
```


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

</details>

This directory contains anything the pipeline downloads, such as assemblies, reads, and databases.


### Aligned genes (`mafft`)

<details markdown="1">
<summary>Output files</summary>

- `aligned_genes/`
  - `busco_genes/`
    - `<gene ID>_aligned.fas`: FASTA files of aligned genes.
  - `core_genes/`
    - `<gene ID>_aligned.fas`: FASTA files of aligned genes.

</details>

FASTA files for each gene extracted from assemblies and aligned for each sample and reference.


### BBMap Sendsketch results

<details markdown="1">
<summary>Output files</summary>

- `sendsketch/`
  - `<Sample ID>.txt`: Table returned by BBmap `sendsketch` with initial identifications of samples.

</details>

Tables with information used to make initial identifications of samples from the BBMap `sendsketch` tool.


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


### Percentage of conserved proteins (POCP)

<details markdown="1">
<summary>Output files</summary>

- `pocp`
  - `<Report ID>_pocp.tsv`: A pairwise matrix of the POCP between all samples and references

</details>

POCP is calculated as a metric to comapare samples to eachother and to references in regards to shared gene content


### BUSCO

<details markdown="1">
<summary>output files</summary>

- `busco/`
  - `short_summary.specific.<busco_db>.<species_name>.fasta.txt`: completeness report in tsv format
  - `<species_name>-<busco_db>-busco.batch_summary.txt`: summarized completeness report in tsv format
  - `<sample id>-<database lineage>-busco`: directory with other busco results

</details>

BUSCO is used to extract gene for phylogenetic analysis of eukaryotes and asses assembly completeness.


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

Various trees are produced by the pipeline to compare the samples to references and to eachother.
To put samples in context of reference genomes and provide data that can be useful in identification, core genes from prokaryotes and BUSCO genes from eukaryotes are used to produce trees with `iqtree2`.
SNPs identified by variant calling are also used to create a tree with `iqtree2` for high-resolution sample comparison.


### Quality control reports

<details markdown="1">
<summary>output files</summary>

- `quality_control/`
  - `multiqc/`
    - `<Report ID>_multiqc>`: MultiQC outputs for samples in each report
  - `nanoplot/`: Nanoplot output reports and plots
  - `quast/`
    - `<Sample ID>`: Quast reports and associated data

</details>

Various tools are used to check reads and assemblies for quality.
The outputs of these tools are compiled using MultiQC.


### Main reports

<details markdown="1">
<summary>output files</summary>

- `reports/`
  - `<Report ID>_report.html`: The primary output report of the pipeline

</details>

This is the primary output of the pipeline, containing the report meant to be understandable by non-bioinformations.


### Main report inputs

<details markdown="1">
<summary>output files</summary>

- `main_report_input/`
  - `<Report ID>_inputs`: A folder containing formatted outputs from the pipeline used in the main report.

</details>

This is the directory used to create the main report for each report group.
It contains selected and renamed outputs from the pipeline present in other output folders.


### Pirate

<details markdown="1">
<summary>output files</summary>

- `pirate/`
  - `<Report ID>_results`: Pirate output

</details>

Pirate is used to identify orthologous gene clusters, which is are used later in the pipeline to create phylogenies of prokaryotes with the maximum number of shared genes without relying on annotations.


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
