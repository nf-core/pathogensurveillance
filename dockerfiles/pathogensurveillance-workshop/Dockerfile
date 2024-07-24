################## BASE IMAGE ######################

FROM nfcore/gitpod:latest

################## METADATA #######################

LABEL about.summary="Docker image for running the pathogensurveillance pipeline in workshops."
LABEL about.author="Zachary S.L. Foster, Camilo H. Parada Rojas"
MAINTAINER Zachary S.L. Foster <zacharyfoster1989@gmail.com>

################## BUILD ######################

USER root

RUN apt-get update && apt-get install ffmpeg libsm6 libxext6 -y

RUN <<EOF
conda install --yes nanoplot=1.41.6 quast=5.2.0 multiqc=1.18 bandage=0.8.1
EOF

COPY reads /data/reads
COPY sample_data.csv /data/sample_data.csv
