# Meta pipeline for HiFi assembly and QC
This repository contains a Snakemake pipeline for the assembly, scaffolding, analysis, annotation and quality control of HiFi-based assemblies.
Although developed for a project in lettuce, the pipeline is designed to work with any organism.
The pipeline will work with both HiFi and ONT data, although the former is required.

## Downloading pipeline
The pipeline can be obtained via:
```bash
git clone https://github.com/dirkjanvw/meta_pipeline_hifi.git
cd meta_pipeline_hifi/
```

### (Optional) updating pipeline
Should you notice that bugs have been fixed on GitHub or a new feature implemented in the pipeline, updating the pipeline is as simple as running the following in the `meta_pipeline_hifi/` directory:
```bash
git pull
```

## Installing dependencies
The pipeline will work on any Linux system where `conda`/`mamba` and `singularity`/`apptainer` are installed.
Please run the following to install `snakemake`:
```bash
mamba create -c conda-forge -c bioconda -n snakemake snakemake=8
```
Then activate the environment:
```bash
conda activate snakemake
```

## Databases
Next to Snakemake, `conda`/`mamba` and `singularity`/`apptainer`, this pipeline depends on the existence of a number of databases.

| Database            | Download instructions                                                                                                                            |
|---------------------|--------------------------------------------------------------------------------------------------------------------------------------------------|
| Helixer model       | [Helixer GitHub](https://github.com/weberlab-hhu/Helixer/blob/main/resources/model_list.csv)                                                     |
| GXDB database       | Follow "Download the database" instructions on [FCS GitHub wiki](https://github.com/ncbi/fcs/wiki/FCS-GX) (I only tested the Cloud instructions) |
| Kraken2 nt database | Download `nt` from [this list](https://benlangmead.github.io/aws-indexes/k2)                                                                     |
| OMA database        | Download `LUCA.h5` from [this list](https://omabrowser.org/oma/current/)                                                                         |

## Configuration
First, the configuration file `config/config.yaml` must be created.
Please use the provided `config/example.yaml` as a template and carefully follow the instructions.

## Running the pipeline

### Available modules
Several modules are available in this pipeline (will be referred to later as `${MODULE}`):
- `assemble`: This module will only assembly the reads into contigs.
- `scaffold`: This module will scaffold the contigs using `ntJoin` against a provided reference genome.
- `analyse`: This module will analyse the assembly for provided genes, sequences and contamination.
- `annotate`: This module will generate a quick-and-dirty annotation of the assembly using `liftoff` and `helixer`.
- `qc`: This module will perform quality control of the scaffolded assembly and the quick-and-dirty annotation.
- `all`: This module will run all the above modules.

### Important parameters
Several important `snakemake` parameters are important when running this pipeline:

| Parameter              | Optionality | Description                                                               |
|------------------------|-------------|---------------------------------------------------------------------------|
| `-n`                   | Optional    | Do a dry-run of the pipeline.                                             |
| `-p`                   | Optional    | Print the shell commands that are being executed.                         |
| `-c`/`--cores`         | Required    | Number of CPUs to use; will be referred to as `${CPU}`                    |
| `--use-conda`          | Required    | Use `conda`/`mamba` to manage dependencies.                               |
| `--conda-prefix`       | Optional    | Path where the `conda` environments will be stored.                       |
| `--use-singularity`    | Required    | Use `singularity`/`apptainer` to manage containers.                       |
| `--singularity-prefix` | Optional    | Path where the `singularity` images will be stored.                       |
| `--resources`          | Optional    | Information about system resources; see below at [Resources](#resources). |

### Resources
The following resources (apart from CPUs) might be heavily used by the pipeline:
- `gbmem`: The amount of memory in GB that RAM-heavy jobs in the pipeline can use; will be referred to as `${MEM}`.
- `helixer`: The number of Helixer jobs in the pipeline can run to at the same time; will be referred to as `${HELIXER}`.
- `pantools`: The number of PanTools jobs in the pipeline can run to at the same time; will be referred to as `${PANTOOLS}`.

### Running the pipeline
As first step, it is always good to do a dry-run to check if everything is set up correctly:
```bash
snakemake ${MODULE} -npc${CPU} --use-conda --use-singularity --resources gbmem=${MEM} helixer=${HELIXER} pantools=${PANTOOLS}
```

If everything is alright, the pipeline can be run:
```bash
snakemake ${MODULE} -pc${CPU} --use-conda --use-singularity --resources gbmem=${MEM} helixer=${HELIXER} pantools=${PANTOOLS}
```

### Reporting
The pipeline can generate an HTML `report.html` file with the most important results:
```bash
snakemake ${MODULE} -c1 --report report.html
```

## Output
Next to the [report](#reporting) generated by Snakemake, the most important outputs of the pipeline are the genome assembly and annotation.
These can be found at the following location (where `${assembly}` is a wildcard for each separate assembly):

```bash
results/${assembly}/2.scaffolding/02.renaming/${assembly}.fa
results/${assembly}/4.annotation/03.combined/${assembly}.gff
```
