#!/bin/bash

set -exuo pipefail

# Check if the test data directory exists
if [ ! -d "test_data" ]; then
    mkdir test_data
fi
cd test_data/

# Download kraken2 nt_core
if [ ! -d "nt" ]; then
    mkdir nt
    cd nt
    wget "https://genome-idx.s3.amazonaws.com/kraken/k2_core_nt_20241228.tar.gz"
    tar xzvf k2_core_nt_20241228.tar.gz
    cd ..
fi

# Download OMA LUCA
if [ ! -f "LUCA.h5" ]; then
    wget "https://omabrowser.org/All/LUCA.h5"
fi

# Download the GXDB
if [ ! -d "gxdb" ]; then
    mkdir gxdb
    cd gxdb
    curl -LO https://github.com/peak/s5cmd/releases/download/v2.0.0/s5cmd_2.0.0_Linux-64bit.tar.gz
    tar -xvf s5cmd_2.0.0_Linux-64bit.tar.gz
    ./s5cmd --no-sign-request cp --part-size 50 --concurrency 50 s3://ncbi-fcs-gx/gxdb/latest/all.* ./
    cd ..
fi

# Return to the main directory
cd ..
