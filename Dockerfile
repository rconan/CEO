FROM nvidia/cuda:9.1-cudnn7-devel
FROM continuumio/anaconda3

RUN apt-get update && apt-get -y install noweb

RUN apt-get -y install libcurl4-openssl-dev

RUN apt-get -y install git

COPY . /tmp/src

RUN cd /tmp/src && git checkout devel_cuda9.1_python3.6 && make all cython
