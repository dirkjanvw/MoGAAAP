# Meta pipeline for HiFi assembly and QC
This repository contains a Snakemake pipeline for the assembly, scaffolding, analysis, annotation and quality control of HiFi-based assemblies.
Although developed for a project in lettuce, the pipeline is designed to work with any organism.
The pipeline will work with both HiFi and ONT data, although the former is required.

## Index
- [Downloading the pipeline](#downloading-pipeline)
- [Installing dependencies](#installing-dependencies)
- [Required databases](#databases)
- [Configuring the pipeline](#configuration)
- [Running the pipeline](#running-the-pipeline)
- [Output](#output)
- [Explaining the pipeline](#explaining-the-pipeline)
- [FAQ](#faq)

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

### Conda/mamba
If not installed already, `conda`/`mamba` can be installed by following these instructions:
```bash
# install miniforge
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
# follow instructions and let it run `conda init`

# set default channels
source ~/.bashrc
conda config --set auto_activate_base false
source ~/.bashrc
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

### Singularity/Apptainer
This installation process requires root access, but typically a server admin can install it for you if it is not already installed.
If not installed already, Singularity/Apptainer can be installed by following the instructions on their website: [apptainer.org](https://apptainer.org/docs/user/latest/quick_start.html).
Make sure to install either Singularity version 3.7 or higher or Apptainer version 1.0 or higher.

**NB**: For running the pipeline, no root access or special permissions are required.
However, the pipeline needs some SIF files that are not included in the repository.
These need to be built using the provided DEF files, which requires `sudo` permissions.
Building these SIF files can be done locally on a personal laptop (or on a server if you happen to have `sudo` permissions).
Navigate to the `workflow/singularity` directory and run the following command for each DEF file (file ending in `.def`):
```bash
sudo singularity build ${NAME}.sif ${NAME}.def
```

### Snakemake
Snakemake can be installed using `conda`/`mamba`:
```bash
mamba create -c conda-forge -c bioconda -n snakemake snakemake=8
```

Then activate the environment before running the pipeline:
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
All configuration of the pipeline is done in the `config/config.yaml` file, and samples are registered in the `config/samples.tsv` file.
All fields to fill in are well-documented in the provided `config/config.yaml` file and should be self-explanatory.

The `config/samples.tsv` has the following columns to fill in (one row per sample):

| Column name     | Description                                                                                                                                                  |
|-----------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `accessionId`   | The accession ID of the sample. This name has to be unique.                                                                                                  |
| `hifi`          | The path to the HiFi reads in FASTQ or FASTA format.                                                                                                         |
| `ont`           | The path to the ONT reads in FASTQ or FASTA format.                                                                                                          |
| `illumina_1`    | The path to the forward Illumina reads in FASTQ format.                                                                                                      |
| `illumina_2`    | The path to the reverse Illumina reads in FASTQ format.                                                                                                      |
| `speciesName`   | A name for the species that is used by Helixer to name novel genes.                                                                                          |
| `taxId`         | The NCBI taxonomy ID of the species.                                                                                                                         |
| `referenceId`   | A unique identifier for the reference genome for which genome (FASTA), annotation (GFF3) and chromosome names are provided in the `config/config.yaml` file. |

Both `config/config.yaml` and `config/samples.tsv` files validated against a built-in schema that throws an error if the files are not correctly filled in.

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
For more information about these modules, see [Explaining the pipeline](#explaining-the-pipeline).

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

## Explaining the pipeline
Assembling a genome from raw data to a final usable resource is a process that is hard to automate.
We believe that this process always necessesitates human curation.
However, large parts can easily be automated, which is why we created this pipeline.
This pipeline performs the assembly, scaffolding and renaming of genomic data as well as an initial quick-and-dirty structural annotation.
Importantly, both genome assembly and annotation are subjected to quality control, providing a direct starting point for the curation of the assembly.
Furthermore, each part of the process (module) can be run separately after which its output can be inspected before continuing to the next step.

```mermaid
graph TD;
    assemble-->scaffold;
    scaffold-->analyse;
    scaffold-->annotate;
    scaffold-->qc;
    annotate-->qc;
```

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
Proper structural genome annotation would take too long and is not a problem that is solved for automation yet.
Therefore, we implemented a "quick-and-dirty" annotation in this pipeline by combing the results of `liftoff` and `helixer`.
`helixer` will run on the GPU if it's available, otherwise it will run on CPU (which is known to be a lot slower).
In case of overlap in features between `liftoff` and `helixer`, we take the `liftoff` annotation.

#### Next steps
Although this module generally runs for the longest time, no visual output is produced; only a final GFF3 file.
This GFF3 file, together with the FASTA file from the scaffold module are the only inputs for the final module: QC.

### QC module

#### Overview
This final quality control module is the most important for human curation of the genome.
The quality control steps in this module can be roughly divided into two categories: individual and grouped.
Individual quality control steps include k-mer completeness (`merqury`), k-mer contamination (`kraken2`), NCBI contamination (`fcs-gx`), adapter contamination (`fcs-adaptor`) and read mapping (`bwa-mem2`).
Grouped quality control steps include BUSCO completeness (`busco`), OMA completeness (`omark`), k-mer distances (`kmer-db`), mash distances (`mash`), minimizer collinearity (`ntsynt`), k-mer phylogeny (`SANS`), k-mer pangenome growth (`pangrowth`), gene pangenome growth (`pantools`) and general statistics.
These groups are meant to give a comparative overview of the assembly and annotation.
Any groups can be defined in the configuration file and a genome may occur in multiple groups.

#### Next steps
The report (see [Reporting](#reporting)) produced by this module is the most useful output of the pipeline for human curation.
It contains visual output for each of the quality control steps performed in this module including a description on how to interpret the results.
Importantly, the qc module does not do any filtering of the assembly or annotation, only reporting.
Next steps could include (but are not limited to) removal of contaminants, discovery of sample swaps, subsetting the input data, etc.

## FAQ

### Q: Pipeline crashes at renaming the chromosomes
A: This issue typically arises when the assembly and reference genome are not collinear because of an evolutionary distance that is too large.
In this case, the pipeline is not able to accurately discern which reference chromosome corresponds to which assembly scaffold.
This lack of collinearity should also be visible in the MUMmerplots created by the `assemble` module.
The only solution in this case would be to choose another (more closely related) reference genome and re-run the `assemble` module to check if this new reference genome is collinear with the assembly before continuing with the `scaffold` module.

### Q: Pipeline crashes mid-run
A: This is intended Snakemake behaviour: if any job fails, the pipeline will only finish current running jobs and then exit.
As for the reason of stopping, please check the log file of the job that failed.
The name of the log file will be printed in the terminal output of the pipeline.
If the error is not clear, please open an issue on this GitHub page.

### Q: The pipeline cannot find software X
A: Make sure that all SIF containers are built (see [Singularity/Apptainer](#singularityapptainer)) and that the pipeline is run with both `--use-conda` and `--use-singularity`.
All dependencies are included in either a `conda` environment or a `singularity` container.

### Contact
If the above information does not answer your question or solve your issue, feel free to open an issue on this GitHub page or send me an email over dirk[dash]jan[dot]vanworkum[at]wur[dot]nl.


