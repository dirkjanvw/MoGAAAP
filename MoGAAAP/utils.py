import click
import os
import shutil
import subprocess
from yaml import dump
from math import log


# Define globally to make sure init has the correct paths
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
WORKFLOW_DIR = os.path.join(BASE_DIR, 'workflow')
CONFIG_DIR = os.path.join(BASE_DIR, 'config')


def run_command(cmd, output_file=None):
    click.secho(f'[INFO ] Running the following command: {cmd}', fg='blue')
    try:
        if output_file:
            with open(output_file, 'w') as outfile:
                subprocess.run(cmd, check=True, stdout=outfile)
        else:
            subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        click.secho(f'[ERROR] Command failed with exit code {e.returncode}', fg='red')
        return False
    return True


def show_ascii_art():
    click.secho("""
        ___  ___      _____   ___    ___    ___  ______
        |  \\/  |     |  __ \\ / _ \\  / _ \\  / _ \\ | ___ \\
        | .  . | ___ | |  \\// /_\\ \\/ /_\\ \\/ /_\\ \\| |_/ /
        | |\\/| |/ _ \\| | __ |  _  ||  _  ||  _  ||  __/
        | |  | | (_) | |_\\ \\| | | || | | || | | || |
        \\_|  |_/\\___/ \\____/\\_| |_/\\_| |_/\\_| |_/\\_|
        """,
        fg='blue')


def init_mogaaap(workdir):
    click.secho(f'[INFO ] Initialising a new MoGAAAP pipeline at {workdir}',
        fg='blue')

    # Copy the workflow directory to the working directory if it doesn't exist yet
    if not os.path.exists(os.path.join(workdir, 'workflow')):
        shutil.copytree(WORKFLOW_DIR, os.path.join(workdir, 'workflow'))
    else:
        click.secho('[ERROR] Workflow directory already exists; exiting', fg='red')
        return

    # Copy the config directory to the working directory if it doesn't exist yet
    if not os.path.exists(os.path.join(workdir, 'config')):
        shutil.copytree(CONFIG_DIR, os.path.join(workdir, 'config'))
    else:
        click.secho('[ERROR] Config directory already exists; exiting', fg='red')

    # Print further instructions
    click.secho(f'[INFO ] Initialised a new MoGAAAP pipeline at {workdir}',
        fg='blue')
    click.secho(f'[INFO ] Please configure a config.yaml in {workdir}/config/',
        fg='blue')
    click.secho('        before running the pipeline by hand or using `MoGAAAP configure`',
        fg='blue')
    click.secho(f'[INFO ] See {workdir}/config/example.yaml',
                fg='blue')
    click.secho('        for an example configuration',
            fg='blue')

    return


def download_databases(workdir, databases):
    click.secho(f'[INFO ] Downloading databases for MoGAAAP: {databases}', fg='blue')
    databases = [database.lower() for database in databases]

    # Check if the working directory exists
    if not os.path.exists(workdir):
        os.makedirs(workdir)

    # Download kraken2 if needed
    if "all" in databases or "kraken2" in databases:
        kraken2_dir = os.path.join(workdir, "core_nt")

        if os.path.exists(kraken2_dir):
            click.secho('[WARN ] kraken2 database already exists; skipping download', fg='yellow')
        else:
            click.secho(f'[INFO ] Downloading kraken2 database to {kraken2_dir}', fg='blue')
            os.makedirs(kraken2_dir)

            kraken2_download_cmd = ['wget', '"https://genome-idx.s3.amazonaws.com/kraken/k2_core_nt_20241228.tar.gz"',
                '-O', os.path.join(kraken2_dir, 'k2_core_nt_20241228.tar.gz')]
            run_command(kraken2_download_cmd)

            kraken2_untar_cmd = ['tar', 'xzf', os.path.join(kraken2_dir, 'k2_core_nt_20241228.tar.gz'),
                '-C', kraken2_dir]
            run_command(kraken2_untar_cmd)

            click.secho(f'[INFO ] kraken2 database downloaded to {kraken2_dir}', fg='blue')

    # Download OMA if needed
    if "all" in databases or "oma" in databases:
        oma_dir = os.path.join(workdir, "oma")

        if os.path.exists(oma_dir):
            click.secho('[WARN ] OMA database already exists; skipping download', fg='yellow')
        else:
            click.secho(f'[INFO ] Downloading OMA database to {oma_dir}', fg='blue')
            os.makedirs(oma_dir)

            oma_download_cmd = ['wget', '"https://omabrowser.org/All/LUCA.h5"',
                '-O', os.path.join(oma_dir, 'LUCA.h5')]
            run_command(oma_download_cmd)

            click.secho(f'[INFO ] OMA database downloaded to {oma_dir}', fg='blue')

    # Download GXDB if needed
    if "all" in databases or "gxdb" in databases:
        gxdb_dir = os.path.join(workdir, "gxdb")

        if os.path.exists(gxdb_dir):
            click.secho('[WARN ] GXDB database already exists; skipping download', fg='yellow')
        else:
            click.secho(f'[INFO ] Downloading GXDB database to {gxdb_dir}', fg='blue')
            os.makedirs(gxdb_dir)

            gxdb_helper_cmd = ['wget', '"https://github.com/peak/s5cmd/releases/download/v2.0.0/s5cmd_2.0.0_Linux-64bit.tar.gz"',
                '-O', os.path.join(gxdb_dir, 's5cmd_2.0.0_Linux-64bit.tar.gz')]
            run_command(gxdb_helper_cmd)

            gxdb_untar_cmd = ['tar', 'xzf', os.path.join(gxdb_dir, 's5cmd_2.0.0_Linux-64bit.tar.gz'),
                '-C', gxdb_dir]
            run_command(gxdb_untar_cmd)

            gxdb_download_cmd = [os.path.join(gxdb_dir, 's5cmd'), '--no-sign-request', 'cp', '--part-size', '50',
                '--concurrency', '50', 's3://ncbi-fcs-gx/gxdb/latest/all.*', gxdb_dir]
            run_command(gxdb_download_cmd)

            click.secho(f'[INFO ] GXDB database downloaded to {gxdb_dir}', fg='blue')

    return


def configure_mogaaap(workdir, samples, reference_fasta, reference_gff,
    mitochondrion, chloroplast, telomere, odb, helixer_model, gxdb, omadb,
    kraken2db):
    click.secho('[INFO ] Configuring the MoGAAAP pipeline', fg='blue')

    # Check if the workflow directory exists
    if not os.path.exists(os.path.join(workdir, 'workflow')):
        click.secho('[ERROR] Workflow directory does not exist', fg='red')
        click.secho('[ERROR] Did you run `MoGAAAP init`?', fg='red')
        return

    # Check if the configfile already exists and move it if it does
    configfile = os.path.join(workdir, 'config', 'config.yaml')
    if os.path.exists(configfile):
        old_config_counter = 1
        while os.path.exists(os.path.join(workdir, 'config', f'config_{old_config_counter}.yaml')):
            old_config_counter += 1
        old_config_location = os.path.join(workdir, 'config', f'config_{old_config_counter}.yaml')

        click.secho(f'[WARN ] Configuration file already exists; moving it to {old_config_location}',
            fg='yellow')
        os.rename(configfile, old_config_location)

    # Check if the samples contains all required columns
    required_columns = ['accessionId', 'haplotypes', 'speciesName', 'taxId', 'referenceId']
    required_oneof_columns = ['assemblyLocation', 'hifi']
    optional_columns = ['annotationLocation', 'ont', 'illumina_1', 'illumina_2', 'hic_1', 'hic_2']
    accession_names = []
    accession_name_col = -1
    reference_names = []
    reference_name_col = -1
    col_counter = 0
    with open(samples) as f:
        header = f.readline().strip().split('\t')
        for column in required_columns:
            if column not in header:
                click.secho(f'[ERROR] Required column "{column}" not found in {samples}',
                    fg='red')
                return
        required_oneof_column_count = 0
        for column in required_oneof_columns:
            if column in header:
                required_oneof_column_count += 1
        if required_oneof_column_count == 0:
            click.secho(f'[ERROR] One of the required columns {", ".join(required_oneof_columns)} not found in {samples}',
                fg='red')
            return
        for column in header:
            if column not in required_columns and column not in required_oneof_columns and column not in optional_columns:
                click.secho(f'[WARN ] Unknown column "{column}" found in {samples}',
                    fg='yellow')
            if column == 'referenceId':
                reference_name_col = col_counter
            if column == 'accessionId':
                accession_name_col = col_counter
            col_counter += 1
        for line in f:
            reference_name = line.strip().split('\t')[reference_name_col]
            if reference_name not in reference_names:
                reference_names.append(reference_name)
            accession_name = line.strip().split('\t')[accession_name_col]
            if accession_name not in accession_names:
                accession_names.append(accession_name)

    # Obtain the name of the reference genome (from the reference FASTA file)
    # and check if it's part of the referenceId column in the samples
    reference_name = os.path.basename(reference_fasta).split('.')[0]
    if reference_name not in reference_names:
        click.secho(f'[WARN ] Reference genome "{reference_name}" not found in {samples}; only found {", ".join(reference_names)}',
            fg='yellow')
    if len(reference_names) > 1:
        click.secho(f'[WARN ] Multiple reference genomes found in {samples}: {", ".join(reference_names)}',
            fg='yellow')
        click.secho(f'[WARN ] Please manually add them to the {configfile} file like "{reference_name}"',
            fg='yellow')

    # Calculate the optimal k-mer size for the genome from the reference length
    collision_rate = 0.001
    reference_length = 0
    chr_names = {}
    chr_num = 1
    with open(reference_fasta) as f:
        for line in f:
            if line.startswith('>'):
                chr_names[chr_num] = line.strip().split()[0][1:]
                chr_num += 1
            else:
                reference_length += len(line.strip())
    kmer_size = int(log(reference_length*(1-collision_rate)/collision_rate)/log(4))

    # Create a telomere motif file containing 100x the telomere motif
    telomerefile = os.path.join(workdir, 'config', 'telomere.fa')
    with open(telomerefile, 'w') as f:
        f.write(f'>telomere\n{telomere*100}\n')

    # Generate the configuration
    configuration = {}
    configuration["samples"] = samples
    configuration["assembler"] = "hifiasm"
    configuration["scaffolder"] = "ntjoin"
    configuration["min_contig_len"] = 10000
    configuration["ntjoin_k"] = 52
    configuration["ntjoin_w"] = 16000
    configuration["reference_genomes"] = {
        reference_name: {
            "genome": reference_fasta,
            "annotation": reference_gff,
            "chromosomes": chr_names
        }}
    configuration["nucl_queries"] = {"telomere": telomerefile}
    configuration["organellar"] = {"mitochondrion": mitochondrion}
    if chloroplast:
        configuration["organellar"]["chloroplast"] = chloroplast
    configuration["telomere_motif"] = telomere
    configuration["helixer_model"] = helixer_model
    configuration["helixer_max_gene_length"] = 64152
    configuration["k_qa"] = kmer_size
    configuration["gxdb"] = gxdb
    configuration["pantools_grouping"] = 3
    configuration["odb"] = odb
    configuration["kraken2_nt"] = kraken2db
    configuration["OMAdb"] = omadb
    configuration["set"] = {"all": accession_names}
    configuration["jvm"] = "-Xmx100g -Xms100g"
    configuration["tmpdir"] = "/dev/shm/"
    configuration["nucmer_maxgap"] = 1000
    configuration["nucmer_minmatch"] = 1000
    configuration["custom_singularity"] = True

    # Write the configuration file
    dump(configuration, open(configfile, 'w'))

    # Print further instructions
    click.secho(f'[INFO ] Configured the MoGAAAP pipeline at {workdir}',
        fg='blue')
    click.secho(f'[INFO ] Please carefully compare the {configfile}',
        fg='blue')
    click.secho(f'        file with the example in {os.path.join(workdir, "config", "example.yaml")}',
        fg='blue')
    click.secho('        for more information on the configuration options',
        fg='blue')
    click.secho('[INFO ] After finishing the configuration, run the pipeline using `MoGAAAP run`',
        fg='blue')


def run_mogaaap(workdir, configfile, reportfile, cores, memory, dryrun, other,
    targets):
    click.secho(f'[INFO ] Running the MoGAAAP pipeline at {workdir}', fg='blue')

    # Check if the configfile exists
    if not os.path.exists(configfile):
        if not os.path.exists(os.path.join(workdir, 'config', 'config.yaml')):
            click.secho('[ERROR] Configuration file config/config.yaml does not exist',
                fg='red')
            click.secho(f'[ERROR] Checked {os.path.abspath(configfile)} and {os.path.abspath(os.path.join(workdir, "config", "config.yaml"))}',
                fg='red')
            click.secho('[ERROR] Did you run `MoGAAAP init`?',
                fg='red')
            return
        else:
            configfile = os.path.abspath(os.path.join(workdir, 'config', 'config.yaml'))

    # Check if Snakemake and Singularity (or Apptainer) are available
    if not shutil.which('snakemake'):
        click.secho('[ERROR] Snakemake is not available', fg='red')
        return
    if not shutil.which('singularity') and not shutil.which('apptainer'):
        click.secho('[ERROR] Singularity or Apptainer is not available',
            fg='red')
        return

    # Build the Snakemake command
    snakemake_cmd = ['snakemake', '--directory', workdir, '--snakefile', os.path.join(workdir, 'workflow', 'Snakefile')]
    if targets:
        snakemake_cmd.extend(targets)
    snakemake_cmd.extend(['--nolock'])
    snakemake_cmd.extend(['--configfile', os.path.abspath(configfile)])
    snakemake_cmd.extend(['--cores', str(cores)])
    snakemake_cmd.extend(['--resources', f'gbmem={memory}'])
    if other:
        snakemake_cmd.extend(other.split(" "))
    if dryrun:
        snakemake_cmd.append('--dryrun')

    # Run the Snakemake command
    if not run_command(snakemake_cmd):
        return

    # If it is a dryrun, return now
    if dryrun:
        return

    # Build the Snakemake command with the report flag
    report_cmd = snakemake_cmd + ['--report', reportfile + '.tmp.html']
    fix_report_cmd = ['sed', '-E', 's/([^l]) h-screen/\1/g',
        reportfile + '.tmp.html']
    remove_tmp_cmd = ['rm', reportfile + '.tmp.html']

    # Run the Snakemake command with the report flag
    if not run_command(report_cmd):
        return
    if not run_command(fix_report_cmd, output_file=reportfile):
        return
    if not run_command(remove_tmp_cmd):
        return

    # Print further instructions
    click.secho('[INFO ] Finished running the MoGAAAP pipeline', fg='blue')
    click.secho(f'[INFO ] Report available at {reportfile}', fg='blue')
    click.secho(f'[INFO ] Final output available at {os.path.join(workdir, "final_output")}',
        fg='blue')
