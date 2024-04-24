# Meta pipeline for HiFi assembly and QC
This repository contains a Snakemake pipeline for the assembly, scaffolding, analysis, annotation and quality control of HiFi-based assemblies.
Although developed for a project in lettuce, the pipeline is designed to work with any organism.
The pipeline will work with both HiFi and ONT data, although the former is required.

## Index
- [Downloading the pipeline](downloading-pipeline)
- [Installing dependencies](installing-dependencies)
- [Required databases](databases)
- [Configuring the pipeline](configuration)
- [Running the pipeline](running-the-pipeline)
- [Output](output)
- [Breaking down the pipeline](breaking-down-the-pipeline)
- [Contact](questions-or-issues)

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
Please use the provided `config/example.yaml` as a template and carefully follow the instructions in the file.

## Running the pipeline

### Available modules
Several modules are available in this pipeline (will be referred to later as `${MODULE}`):
- `assemble`: This module will only assembly the reads into contigs.
- `scaffold`: This module will scaffold the contigs using `ntJoin` against a provided reference genome.
- `analyse`: This module will analyse the assembly for provided genes, sequences and contamination.
- `annotate`: This module will generate a quick-and-dirty annotation of the assembly using `liftoff` and `helixer`.
- `qc`: This module will perform quality control of the scaffolded assembly and the quick-and-dirty annotation.
- `all`: This module will run all the above modules.

It is advisable to run the pipeline module by module for a new set of assemblies and critically look at the results of each module before continuing.
All modules except for `annotate` have visual output that can be inspected in an HTML report file (see at [Reporting](#reporting)).

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

## Breaking down the pipeline

### Assemble module

#### Overview
In the assemble module, HiFi reads are assembled using `hifiasm`.
If ONT reads are given, these are used in the assembly process using the `--ul` parameter of `hifiasm`.
Since the output of `hifiasm` is a GFA file, we next convert this to a FASTA file.
Finally, we produce an alignment of the assembly against the provided reference genome using `nucmer`.
To prevent spurious alignments, we slightly increased the `-l` and `-g` parameter of `nucmer`.

#### Next steps
Since the next step after assembly is the scaffolding process, there has to be collinearity between the assembly and the reference genome.
This can be checked in the dotplot created from the `nucmer` alignment.
If there is no sign of collinearity between the two, reference-guided scaffolding will be impossible.
The only solution in that case would be to choose another (more closely related) reference genome.

### Scaffold module

#### Overview
Scaffolding is performed using `ntJoin`, which uses a minimizer-based reference-guided scaffolding method.
If by visual inspection collinearity between the assembly and reference genome was found, the scaffolding module generally runs without issues.
Should any error occur, please read the corresponding log file of the step that produced the error.
In most cases, the error may be resolved by choosing different values for the `ntjoin_k` and `ntjoin_w` in the config file.

After scaffolding, the sequences in the scaffolded assembly are renamed to reflect their actual chromosome names according to the reference genome.
Finally, `nucmer` is run again to produce an alignment plot for visual inspection of the scaffolding process.

#### Next steps
As the assembly as outputted by this module is used as starting point for the analyse, annotate and qc modules, it is crucial it matches the expectations in terms of size and chromosome number.
Please carefully look at the `nucmer` alignment plot to check that the assembly looks as expected before continuing to a next module.

### Analyse module

#### Overview
The purpose of the analyse module is to create a genome-wide overview of the newly created assembly.
It does this by running several user-defined queries against the assembly using BLAST and creating a `circos` plot.

#### Next steps
Circos plots are notoriously hard to automate and that is no different for this pipeline.
Although a `circos` plot should always be produced, it typically doesn't look quite right yet.
Feel free to copy the config files produced by this pipeline and adjust to your own plotting needs.

### Annotate module

#### Overview
Asdf

#### Next steps
Asdf

### QC module

#### Overview
Asdf

#### Next steps
Asdf

## Questions or issues
In case of any questions or issues with the pipeline, feel free to open an issue on this GitHub page or send me an email over dirk[dash]jan[dot]vanworkum[at]wur[dot]nl
