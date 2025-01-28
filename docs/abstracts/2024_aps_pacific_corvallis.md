# Pathogensurveillance: An automated computational pipeline for identification, population genomics, and monitoring of pathogens

ZSL Foster1, M Sudermann2, C Parada-Rojas2, F Iruegas-Bocardo2, R Alcala-Briseno2, JH Chang2, and NJ Grünwald1

1 Horticultural Crops Disease and Pest Management Research Unit, USDA ARS, Corvallis, Oregon, USA

2 Dept. Botany and Plant Pathology, Oregon State University, Corvallis, Oregon, USA

Rapid and automated analysis of plant pathogen genome sequences is essential for more effective responses to disease outbreaks. To account for the computational complexities of whole genome sequence analysis, we developed the pathogensurveillance pipeline to identify individuals or groups of pathogens and their diversity. Pathogensurveillance analyzes raw reads derived from eukaryotic or prokaryotic pathogens. The pipeline automatically determines the taxonomic placement, finds the closest reference genome sequence, maps sequence reads, calls variants, identifies core genes, and reports back phylogenetic relationships and identifications. References are selected automatically from the NCBI RefSeq database or can be supplied by the user. For each user-defined group of samples, pathogensurveillance produces an HTML report with summary statistics, interactive figures, and a static PDF report. Additionally, users can access the intermediate files in order to conduct further analyses. This pipeline can be executed on any Linux computer or cloud computing environment. Our pipeline automates and accelerates analysis of whole genome sequence data, which is essential for rapid responses to disease outbreaks, and also allows genomic analysis by users with less experience in biocomputing.
