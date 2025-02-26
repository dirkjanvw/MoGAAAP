import click
import os
import shutil


BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
WORKFLOW_DIR = os.path.join(BASE_DIR, 'workflow')
CONFIG_DIR = os.path.join(BASE_DIR, 'config')


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
    click.secho(f'[INFO ] Initialising a new MoGAAAP pipeline at {workdir}', fg='blue')

    # Copy the workflow directory to the working directory if it doesn't exist yet
    if not os.path.exists(os.path.join(workdir, 'workflow')):
        shutil.copytree(WORKFLOW_DIR, os.path.join(workdir, 'workflow'))
    else:
        click.secho('[WARN ] Workflow directory already exists', fg='yellow')

    # Copy the config directory to the working directory if it doesn't exist yet
    if not os.path.exists(os.path.join(workdir, 'config')):
        shutil.copytree(CONFIG_DIR, os.path.join(workdir, 'config'))
    else:
        click.secho('[WARN ] Config directory already exists', fg='yellow')

    # Print further instructions
    click.secho(f'[INFO ] Initialised a new MoGAAAP pipeline at {workdir}', fg='blue')
    click.secho(f'[INFO ] Please configure the YAML and TSV files in {workdir}/config/ before', fg='blue')
    click.secho('        running the pipeline by hand or using `mogaaap configure`', fg='blue')

    return


def configure_mogaaap(workdir, configfile):
    click.secho('Configuring the MoGAAAP pipeline', fg='blue')


def run_mogaaap(workdir, configfile, reportfile, cores, memory, dryrun, other, targets):
    click.secho('Running the MoGAAAP pipeline', fg='blue')

    # Check if the workflow directory exists
    if not os.path.exists(os.path.join(workdir, 'workflow')):
        click.secho('[ERROR] Workflow directory does not exist', fg='red')
        return

    # Check if Snakemake and Singularity (or Apptainer) are available
    if not shutil.which('snakemake'):
        click.secho('[ERROR] Snakemake is not available', fg='red')
        return
    if not shutil.which('singularity') and not shutil.which('apptainer'):
        click.secho('[ERROR] Singularity or Apptainer is not available', fg='red')
        return

    # Build the Snakemake command
    snakemake_cmd = ['snakemake', '--snakefile', os.path.join(workdir, 'workflow', 'Snakefile')]
    snakemake_cmd.extend(['--configfile', configfile])
    snakemake_cmd.extend(['--cores', str(cores)])
    snakemake_cmd.extend(['--resources', f'gbmem={memory}'])
    if other:
        snakemake_cmd.extend(other)
    if dryrun:
        snakemake_cmd.append('--dryrun')
    if targets:
        snakemake_cmd.extend(targets)

    # Run the Snakemake command
    click.secho(f'[INFO ] Running the following command: {" ".join(snakemake_cmd)}', fg='blue')
