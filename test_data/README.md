# Test data for MoGAAAP
This directory contains test data for the MoGAAAP project.
It can be used to test the functionality of the MoGAAAP project.
Please follow the instructions below to run this test.

## Obtain the pipeline
To obtain the pipeline, please clone the MoGAAAP repository from GitHub.
```bash
git clone https://github.com/dirkjanvw/MoGAAAP.git
cd MoGAAAP
```

## Download the test data
Run the following command to download all test data into the `test_data` directory.

```bash
bash test_data/download_test_data.sh
```

## Download the databases
Run the following command to download all databases into the `test_data` directory.
Normally, it would not be recommended to download them into the `test_data` directory due to speed, but for the purpose of this test, it is fine.

> [!WARNING]
> This script will download >900GB of data.
> Make sure you have enough disk space available.

```bash
bash test_data/download_databases.sh
```

## Run the pipeline
To run the pipeline, execute the following command.
```bash
snakemake --configfile test_data/config.yaml
```
