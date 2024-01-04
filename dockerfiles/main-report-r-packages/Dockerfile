################## BASE IMAGE ######################

FROM condaforge/mambaforge:23.1.0-4

################## METADATA #######################

LABEL base_image="mambaforge"
LABEL version="23.1.0-4"
LABEL about.summary="A test environment to generate the output report for pathogensurveillance."
LABEL software="main-report-r-packages"

LABEL about.author="Zachary S.L. Foster"

################## MAINTAINER ######################
MAINTAINER Zachary S.L. Foster <zacharyfoster1989@gmail.com>

USER root

RUN conda install -y conda-forge::quarto conda-forge::r-knitr

