FROM continuumio/miniconda3:4.7.12
MAINTAINER Sangram Keshari Sahu <sangramsahu15@gmail.com>

COPY envs/main.yml .
RUN \
   apt-get install procps ttf-dejavu -y \
   && conda env update -n root -f main.yml \
   && conda clean -a
