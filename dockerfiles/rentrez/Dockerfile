################## BASE IMAGE ######################

FROM eddelbuettel/r2u:22.04

################## METADATA #######################

LABEL about.summary="A container for the rentrez package used to parse input data in the pathogensurvellance pipeline."
LABEL about.author="Zachary S.L. Foster"
MAINTAINER Zachary S.L. Foster <zacharyfoster1989@gmail.com>

################## BUILD ######################

RUN <<EOF
install.r \
    rentrez
EOF
