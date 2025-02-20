#!/bin/bash

set -exuo pipefail

# Check if the test data directory exists
if [ ! -d "test_data" ]; then
    mkdir test_data
fi
cd test_data/

# Download the test data
wget "ftp://download.big.ac.cn/gsa2/CRA008584/CRR591691/CRR591691.fastq.gz" #43_elk_1
wget "ftp://download.big.ac.cn/gsa2/CRA008584/CRR591692/CRR591692.fastq.gz" #44_ket_10
wget "ftp://download.big.ac.cn/gsa2/CRA008584/CRR591693/CRR591693.fastq.gz" #45_meh_0
wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz" #Arabidopsis thaliana Col-0 genome assembly
wget -O Araport11_GFF3_genes_transposons.current.gff.gz "https://www.arabidopsis.org/api/download-files/download?filePath=Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.current.gff.gz" #Arabidopsis thaliana Col-0 genome annotation
wget -O mitochondrion.fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&id=NC_037304.1&report=fasta&format=text" #Arabidopsis thaliana Col-0 mitochondrion genome
wget -O chloroplast.fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&id=NC_000932.1&report=fasta&format=text" #Arabidopsis thaliana Col-0 chloroplast genome

# Unzip the reference files
gunzip GCF_000001735.4_TAIR10.1_genomic.gff.gz Araport11_GFF3_genes_transposons.current.gff.gz

# Rename the sequence identifiers of the GFF file
sed 's/^Chr1/NC_003070.9/g' Araport11_GFF3_genes_transposons.current.gff |\
    sed 's/^Chr2/NC_003071.7/g' |\
    sed 's/^Chr3/NC_003074.8/g' |\
    sed 's/^Chr4/NC_003075.7/g' |\
    sed 's/^Chr5/NC_003076.8/g' |\
    sed 's/^ChrM/NC_037304.1/g' |\
    sed 's/^ChrC/NC_000932.1/g' > Araport11_GFF3_genes_transposons.current.renamed.gff

# Return to the main directory
cd ..
