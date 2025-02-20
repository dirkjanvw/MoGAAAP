#!/bin/bash

set -exuo pipefail

# Check if the test data directory exists
if [ ! -d "test_data" ]; then
    mkdir test_data
fi
cd test_data/

# Download Helixer model
wget -O land_plant_v0.3_a_0080.h5 "https://zenodo.org/records/10836346/files/land_plant_v0.3_a_0080.h5?download=1"

# Download kraken2 nt_core
wget "https://genome-idx.s3.amazonaws.com/kraken/k2_core_nt_20241228.tar.gz"

# Download OMA LUCA
wget "https://omabrowser.org/All/LUCA.h5"

# Download the GXDB
curl -LO https://github.com/peak/s5cmd/releases/download/v2.0.0/s5cmd_2.0.0_Linux-64bit.tar.gz
tar -xvf s5cmd_2.0.0_Linux-64bit.tar.gz
./s5cmd --no-sign-request cp --part-size 50 --concurrency 50 s3://ncbi-fcs-gx/gxdb/latest/all.* ./

# Return to the main directory
cd ..
