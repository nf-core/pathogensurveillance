################## BASE IMAGE ######################

FROM nfcore/gitpod:latest

################## METADATA #######################

LABEL about.summary="Docker image for running the pathogensurveillance pipeline in workshops."
LABEL about.author="Zachary S.L. Foster, Camilo H. Parada Rojas"
MAINTAINER Zachary S.L. Foster <zacharyfoster1989@gmail.com>

################## BUILD ######################

RUN sudo apt-get update && sudo apt-get install ffmpeg libsm6 libxext6 -y

RUN <<EOF
conda create -n qc --yes bioconda::nanoplot bioconda::quast bioconda::multiqc bioconda::bandage nf-core
EOF

COPY reads /data/reads
COPY sample_data.csv /data/sample_data.csv
