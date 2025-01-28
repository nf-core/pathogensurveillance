FROM ubuntu
MAINTAINER Zachary Foster <zacharyfoster1989@@gmail.com>

USER root
RUN apt-get update && apt-get install -y mash khmer && apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/*
