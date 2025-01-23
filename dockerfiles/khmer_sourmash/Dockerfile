################## BASE IMAGE ######################

FROM condaforge/mambaforge:23.1.0-4

################## METADATA #######################

LABEL base_image="mambaforge"
LABEL about.summary="A combination of the tools bbmap and khmer."
LABEL software="bbmap-khmer"
LABEL about.tags="Genomics"

################## MAINTAINER ######################
MAINTAINER Zachary S.L. Foster <zacharyfoster1989@gmail.com>

USER root

RUN conda install --yes python=3.8.19 bioconda::khmer=3.0.0a3 bioconda::sourmash=4.8.3 --channel conda-forge --channel bioconda --channel defaults --strict-channel-priority
