################## BASE IMAGE ######################

FROM continuumio/miniconda3

################## METADATA #######################

LABEL about.summary="Docker image for running the pathogensurveillance pipeline in workshops."
LABEL about.author="Zachary S.L. Foster, Camilo H. Parada Rojas"
MAINTAINER Zachary S.L. Foster <zacharyfoster1989@gmail.com>

################## BUILD ######################

RUN <<EOF
apt update
apt install --yes less vim
EOF

RUN <<EOF
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create --yes --name nf-core python=3.12 nf-core nextflow
EOF

# Add conda bin to path
ENV PATH /opt/conda/bin:$PATH

# Set up environment
RUN echo "conda activate nf-core" >> ~/.bashrc

# Set entrypoint to bash
# ENTRYPOINT ["bash", "-l", "-c"]

#RUN <<EOF
#nextflow run https://github.com/grunwaldlab/pathogensurveillance -profile conda,aps --outdir test_out -resume -latest
#EOF
