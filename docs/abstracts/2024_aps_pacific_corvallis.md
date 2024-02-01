# Pathogensurveillance: An automated computational pipeline for identification, population genomics, and monitoring of pathogens

Rapid and automated analysis of plant pathogen genomes will allow for more effective responses to emerging pathogens.
We have developed the computational pipeline pathogensurveillance to analyze populations of plant pathogens.
Pathogensurveillance takes raw reads from whole-genome sequencing of unknown organisms, eukaryote or prokaryote, and creates a report with interactive plots and tables.
Currently, pathogensurveillance has two main goals: robust identification of organisms and the analysis of genetic diversity within groups of organisms.
Sendsketch is used for a rapid initial identification, which is used to automatically select and download references from NCBI RefSeq.
Raw reads are then assembled, annotated, and used to construct a core genome phylogeny containing both the unknown inputs and the RefSeq references.
Fine-scale genetic diversity within the input samples is inferred by mapping reads to a reference and calling variants.
The reference can be supplied by the user or selected automatically from RefSeq.
The resulting variants are used to construct a SNP tree and a minimum spanning network.
For each user-defined group of samples, pathogensurveillance produces an HTML report with interactive figures and a static PDF report.
All intermediate files are also available, such as mapped reads, annotated genome assemblies, and tree files.
Pathogensurveillance can be easily installed on and take full advantage of any Linux system, including laptops, desktops, high-performance computing clusters, and commercial cloud providers.
This software provides a way for non-bioinformatitions to rapidly analyze raw reads from unknown organisms and provides researchers with starting point for question-specific analyses.


