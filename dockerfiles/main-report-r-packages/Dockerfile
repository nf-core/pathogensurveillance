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
    texlive-base \
    texlive-xetex \
    lmodern \
    texlive-fonts-recommended \
    texlive-fonts-extra
rm -rf /var/lib/apt/lists/*
EOF

RUN <<EOF
install.r \
    shiny \
    jsonlite \
    htmltools \
    remotes \
    renv \
    rmarkdown \
    quarto \
    dplyr \
    readr \
    knitr \
    readr \
    purrr \
    yaml \
    phylocanvas \
    ape \
    magrittr \
    pheatmap \
    tidyverse \
    palmerpenguins \
    ade4 \
    adegenet \
    poppr \
    ggtree \
    igraph \
    visNetwork \
    phangorn \
    ggplot2 \
    ggnewscale \
    kableExtra \
    plotly \
    webshot2 \
    metacoder \
    ggdendro
EOF

RUN <<EOF
install.r devtools
Rscript -e 'devtools::install_github("grunwaldlab/psminer")'
EOF

RUN <<EOF
curl -LO https://quarto.org/download/latest/quarto-linux-amd64.deb
gdebi --non-interactive quarto-linux-amd64.deb
EOF



