################## BASE IMAGE ######################

FROM nfcore/gitpod:latest

################## METADATA #######################

LABEL about.summary="Docker image for running the pathogensurveillance pipeline in workshops."
LABEL about.author="Zachary S.L. Foster, Camilo H. Parada Rojas"
MAINTAINER Zachary S.L. Foster <zacharyfoster1989@gmail.com>

################## BUILD ######################

USER root

# Install nanoplot, quast, bandage, multiqc

RUN <<EOF
    nextflow run https://github.com/grunwaldlab/pathogensurveillance -profile aps,docker -resume -latest
EOF
