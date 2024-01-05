################## BASE IMAGE ######################

FROM eddelbuettel/r2u:22.04

################## METADATA #######################

LABEL about.summary="A test environment to generate the output report for pathogensurveillance."
LABEL about.author="Zachary S.L. Foster"
MAINTAINER Zachary S.L. Foster <zacharyfoster1989@gmail.com>

################## BUILD ######################

RUN <<EOF
apt-get update && apt install -y --no-install-recommends \
    pandoc \
    pandoc-citeproc \
    curl \
    gdebi-core \
    && rm -rf /var/lib/apt/lists/*
EOF
    
RUN <<EOF 
install.r \
    shiny \
    jsonlite \
    ggplot2 \
    htmltools \
    remotes \
    renv \
    knitr \
    rmarkdown \
    quarto
EOF

RUN curl -LO https://quarto.org/download/latest/quarto-linux-amd64.deb
RUN gdebi --non-interactive quarto-linux-amd64.deb

CMD ["bash"]

