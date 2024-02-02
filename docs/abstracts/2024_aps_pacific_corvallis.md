# Pathogensurveillance: An automated computational pipeline for identification, population genomics, and monitoring of pathogens

ZSL Foster1, M Sudermann2, C Parada-Rojas2, F Iruegas-Bocardo2, R Alcala-Briseno2, JH Chang2, and NJ Grünwald1

1 Horticultural Crops Disease and Pest Management Research Unit, USDA ARS, Corvallis, Oregon, USA

2 Dept. Botany and Plant Pathology, Oregon State University, Corvallis, Oregon, USA

Rapid and automated analysis of plant pathogen genomes will allow for more effective responses to emerging pathogens. However, computational analysis of whole genome sequences is currently complex. We developed the computational pipeline pathogensurveillance to diagnose individuals or communities of pathogens and variants within populations. Pathogensurveillance takes raw reads from whole-genome sequences of unknown eukaryotic or prokaryotic pathogens. The pipeline automatically processes sequence read data and determines taxonomic placement using sendsketch, finds the closest reference genome, maps sequence reads to reference, calls variants and reports back the population genetic and phylogenetic species placement. The reference can be supplied by the user or selected automatically from RefSeq. For each user-defined group of samples, pathogensurveillance produces an HTML report with interactive figures and a static PDF report. All intermediate files are also available for further analysis, including vcf, annotated genome assemblies, and tree files. PathogenSurveillance allows for rapid prototyping, reproducibility, portability, and parallelism. Furthermore, this pipeline can be executed in any linux computer or cloud. Our pipeline promises to accelerate analysis of whole genome sequence data and allow analysis by inexperienced users.



