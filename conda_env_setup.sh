#!/bin/bash
envName=$1
conda create -n $envName python=3.7.4 && \
conda install -c conda-forge -n $envName python-magic && \
conda install -c conda-forge -n $envName scons && \
conda install -c bioconda -n $envName samtools=1.10 bwa
