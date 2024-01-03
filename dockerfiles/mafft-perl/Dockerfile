################## BASE IMAGE ######################

FROM condaforge/mambaforge:23.1.0-4

################## METADATA #######################

LABEL base_image="mambaforge"
LABEL version="23.1.0-4"
LABEL about.summary="A combination of the sequence alignment program MAFFT and Perl."
LABEL software="mafft-perl"
LABEL software.version="7.520"
LABEL about.tags="Genomics"
LABEL about.home="https://mafft.cbrc.jp/alignment/software/"
LABEL about.software="https://mafft.cbrc.jp/alignment/software/"
LABEL about.documentation="https://mafft.cbrc.jp/alignment/software/"
LABEL about.license="BSD"
LABEL about.author="Zachary S.L. Foster"

################## MAINTAINER ######################
MAINTAINER Zachary S.L. Foster <zacharyfoster1989@gmail.com>

USER root

RUN conda install bioconda::mafft=7.520 bioconda::perl-bioperl=1.7.8 conda-forge::parallel=20230522

