> [!WARNING]
> This pipeline is still under development and should not be used for production purposes.
> See the [TODO](#todo) section for more information.

# Meta pipeline for HiFi assembly and QC
This repository contains a Snakemake pipeline for the assembly, scaffolding, analysis, annotation and quality control of HiFi-based assemblies.
Although developed for a project in lettuce, the pipeline is designed to work with any organism.
The pipeline will work with both HiFi and ONT data, although the former is preferred.

## Setup
The pipeline will work on any Linux system where `conda`/`mamba` and `singularity`/`apptainer` are installed.
Please run the following to install `snakemake`:
```bash
mamba create -c conda-forge -c bioconda -n snakemake snakemake=8
```
Then activate the environment:
```bash
conda activate snakemake
```

## Configuration
First, the configuration file `config/config.yaml` must be created.
Please use the provided `config/example.yaml` as a template.

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

## TODO
The following tasks are still to be done before the pipeline is finished:
- [x] Add the annotation module.
- [x] Add the QC module.
- [x] Add the renaming of the scaffolds to their proper chromosome names.
- [x] Add the generation of a report.
- [x] Replace compleasm by BUSCO v5.7.0
- [x] Do a full test of the pipeline.
