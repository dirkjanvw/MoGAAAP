# Test data for MoGAAAP
This directory contains test data for the MoGAAAP project.
It can be used to test the functionality of the MoGAAAP project.
Please follow the instructions below to run this test.

## Obtain the pipeline
To obtain the pipeline, please clone the MoGAAAP repository from GitHub:

```bash
git clone https://github.com/dirkjanvw/MoGAAAP.git
cd MoGAAAP
```

## Download the test data
Run the following command to download all test data into the `test_data` directory:

```bash
bash test_data/download_test_data.sh
```

## Download the databases
Run the following command to download all databases into the `test_data` directory.
Normally, it is recommended to download the databases on a disk with fast I/O capabilities, such as an SSD or `/dev/shm/` (if available).
However, for the purpose of this test, it is fine to download them into the `test_data` directory on a regular disk.

> [!WARNING]
> This script will download >900GB of data.
> Make sure you have enough disk space available.

```bash
bash test_data/download_databases.sh
```

## Run the pipeline
To run the pipeline, execute the following command:
```bash
snakemake -s MoGAAAP/workflow/Snakefile --configfile test_data/config.yaml
```

On a system with 128 CPUs, 1TB of memory and no GPUs, this pipeline will take approximately X hours to complete.

## Check the results
After running the pipeline, you can create the HTML report with the following command:
```bash
snakemake -s MoGAAAP/workflow/Snakefile --configfile test_data/config.yaml --report report.html
```

Additionally, the following output files will be present in the `final_output` directory:
```bash
final_output/
├── 43_elk_1.chloroplast.fa
├── 43_elk_1.contigs.fa
├── 43_elk_1.full.clean.gff
├── 43_elk_1.full.coding.gff
├── 43_elk_1.full.fa
├── 43_elk_1.full.gff
├── 43_elk_1.mitochondrion.fa
├── 43_elk_1.nuclear.fa
├── 44_ket_10.chloroplast.fa
├── 44_ket_10.contigs.fa
├── 44_ket_10.full.clean.gff
├── 44_ket_10.full.coding.gff
├── 44_ket_10.full.fa
├── 44_ket_10.full.gff
├── 44_ket_10.mitochondrion.fa
├── 44_ket_10.nuclear.fa
├── 45_meh_0.chloroplast.fa
├── 45_meh_0.contigs.fa
├── 45_meh_0.full.clean.gff
├── 45_meh_0.full.coding.gff
├── 45_meh_0.full.fa
├── 45_meh_0.full.gff
├── 45_meh_0.mitochondrion.fa
└── 45_meh_0.nuclear.fa
```
