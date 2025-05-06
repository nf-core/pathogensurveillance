# Dockerfile to create container with Cell Ranger v8.0.0 and bcl2fastq v2.20.0
# Push to nfcore/ncbi-datasets-cli:<VER>

FROM condaforge/mambaforge:23.1.0-4
LABEL base_image="mambaforge"
LABEL version="23.1.0-4"
LABEL about.summary="NCBI Datasets is a new resource that lets you easily gather data from across NCBI databases."
LABEL software="ncbi-datasets-cli"
LABEL software.version="15.12.0"
LABEL about.tags="Genomics"
LABEL about.home="https://www.ncbi.nlm.nih.gov/datasets"
LABEL about.software="https://www.ncbi.nlm.nih.gov/datasets"
LABEL about.documentation="https://www.ncbi.nlm.nih.gov/datasets/docs/"
LABEL about.license="public_domain"
LABEL about.licence_file="https://github.com/ncbi/datasets/blob/master/LICENSE.md"
LABEL about.author="NCBI"
LABEL extra.identifiers.biotools="ncbi_datasets"
MAINTAINER Zachary S.L. Foster <zacharyfoster1989@gmail.com>

USER root

RUN conda install -c conda-forge ncbi-datasets-cli=16.0.0 unzip bioconda::samtools=1.18

