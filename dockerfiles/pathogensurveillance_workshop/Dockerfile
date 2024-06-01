################## BASE IMAGE ######################

FROM continuumio/miniconda3

################## METADATA #######################

LABEL about.summary="A test environment to generate the output report for pathogensurveillance."
LABEL about.author="Zachary S.L. Foster"
MAINTAINER Zachary S.L. Foster <zacharyfoster1989@gmail.com>

################## BUILD ######################

RUN <<EOFV
apt update
apt install less vim
EOF

RUN <<EOFV
conda create --yes --name nf-core python=3.12 nf-core nextflow
EOF

RUN <<EOFV
nextflow run https://github.com/grunwaldlab/pathogensurveillance -profile conda --input https://raw.githubusercontent.com/grunwaldlab/pathogensurveillance/master/test/data/metadata/xanthomonas_small.csv --outdir test_out --download_bakta_db -resume -latest
EOF
