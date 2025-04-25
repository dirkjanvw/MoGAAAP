#!/bin/bash
set -o pipefail

# Update taxonomy for krona
cd $CONDA_PREFIX/opt/krona/taxonomy
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
ktUpdateTaxonomy.sh --only-build
