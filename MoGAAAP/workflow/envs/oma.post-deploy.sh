#!/bin/bash
set -o pipefail

# Update NCBI location
find $CONDA_PREFIX -type f -name "ncbiquery.py" -exec sed -i 's/ncbi.nih/ncbi.nlm.nih/g' {} +

# Update ete3 toolkit database
python3 - <<EOF
from ete3 import NCBITaxa
ncbi = NCBITaxa()
ncbi.update_taxonomy_database()
EOF
