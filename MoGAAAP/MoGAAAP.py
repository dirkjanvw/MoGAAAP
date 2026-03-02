#!/usr/bin/env python3
# Written by Dirk-Jan van Workum (2025)
# This is a wrapper around the MoGAAAP Snakemake workflow
#
# Usage: mogaaap [subcommand] [options]
#
# Subcommands:
# - init: Initialise a new MoGAAAP pipeline
# - download_databases: Download the necessary databases for the pipeline if needed
# - configure: Configure the MoGAAAP pipeline
# - run: Run the MoGAAAP pipeline and create a report

import click
import subprocess
from importlib.metadata import version
import multiprocessing, psutil, os
from .utils import show_ascii_art, init_mogaaap, configure_mogaaap, run_mogaaap, download_databases as utils_download_databases


class OrderedGroup(click.Group):
    def __init__(self, name=None, commands=None, **attrs):
        super().__init__(name, commands, **attrs)
        self.order = []

    def command(self, *args, **kwargs):
        def decorator(f):
            cmd = super(OrderedGroup, self).command(*args, **kwargs)(f)
            self.order.append(cmd.name)
            return cmd
        return decorator

    def list_commands(self, ctx):
        return self.order


def get_version():
    """Obtain the git version via git, if available."""
    base = version("MoGAAAP")
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            stderr=subprocess.DEVNULL,
            cwd=os.path.dirname(__file__)
        ).decode().strip()
        return f"{base}, {commit}"
    except Exception:
        return base


@click.group(cls=OrderedGroup, context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(version=get_version(), prog_name="MoGAAAP")
def cli():
    """This is a wrapper script around the MoGAAAP Snakemake workflow."""
    pass


@cli.command('init',
    short_help='Initialise a new MoGAAAP pipeline')
@click.option('--workdir', '-d',
    default='.',
    show_default=True,
    type=click.Path(exists=False),
    help='Working directory for MoGAAAP')
def init(workdir):
    """Initialise a new MoGAAAP pipeline"""

    workdir = os.path.abspath(workdir)

    init_mogaaap(workdir)


def validate_databases(ctx, param, value):
    """
    Validate that the databases are valid (i.e. one of: all, kraken2, OMA, GXDB)
    """

    valid_databases = ['all', 'kraken2', 'oma', 'gxdb']

    if not value:
        value = ['all']

    for database in value:
        if database.lower() not in valid_databases:
            raise click.BadParameter(f'Invalid database: {database}')

    return value


@cli.command('download_databases',
    short_help='Download the necessary databases for MoGAAAP')
@click.option('--workdir', '-d',
    default='databases',
    show_default=True,
    help='Target directory for the databases')
@click.argument('databases',
    nargs=-1,
    required=False,
    callback=validate_databases)
def download_databases(workdir, databases):
    """Download databases for MoGAAAP
    Possible DATABASES are: all, kraken2, OMA, GXDB (default: all).
    """

    workdir = os.path.abspath(workdir)

    utils_download_databases(workdir, databases)


@cli.command('configure',
    short_help='Configure the MoGAAAP pipeline')
@click.option('--workdir', '-d',
    default='.',
    show_default=True,
    type=click.Path(exists=True),
    help='Working directory for MoGAAAP')
@click.option('--samples', '-s',
    required=True,
    type=click.Path(exists=True),
    help='TSV file with sample information')
@click.option('--reference-fasta', '-f',
    required=True,
    type=click.Path(exists=True),
    help='Reference FASTA file')
@click.option('--reference-gff', '-g',
    required=True,
    type=click.Path(exists=True),
    help='Reference GFF file')
@click.option('--mitochondrion', '-m',
    type=click.Path(exists=True),
    help='Mitochondrion FASTA file')
@click.option('--chloroplast', '-c',
    type=click.Path(exists=True),
    help='Chloroplast FASTA file')
@click.option('--telomere', '-t',
    default='CCCTAAA',
    show_default=True,
    help='Telomere motif')
@click.option('--odb', '-b',
    help='ODB name for BUSCO')
@click.option('--helixer-model', '-e',
    help='Helixer model name (land_plant, vertebrate, invertebrate, fungi)')
@click.option('--gxdb', '-x',
    type=click.Path(exists=True),
    help='GXDB database location')
@click.option('--omadb', '-o',
    type=click.Path(exists=True),
    help='OMA database location')
@click.option('--kraken2db', '-k',
    type=click.Path(exists=True),
    help='kraken2 database location')
def configure(workdir, samples, reference_fasta, reference_gff, mitochondrion,
    chloroplast, telomere, odb, helixer_model, gxdb, omadb, kraken2db):
    """Configure the MoGAAAP pipeline"""

    workdir = os.path.abspath(workdir)
    samples = os.path.abspath(samples)
    reference_fasta = os.path.abspath(reference_fasta)
    reference_gff = os.path.abspath(reference_gff)
    if mitochondrion:
        mitochondrion = os.path.abspath(mitochondrion)
    if chloroplast:
        chloroplast = os.path.abspath(chloroplast)
    if gxdb:
        gxdb = os.path.abspath(gxdb)
    if omadb:
        omadb = os.path.abspath(omadb)
    if kraken2db:
        kraken2db = os.path.abspath(kraken2db)

    configure_mogaaap(workdir, samples, reference_fasta, reference_gff,
        mitochondrion, chloroplast, telomere, odb, helixer_model, gxdb, omadb,
        kraken2db)


def validate_targets(ctx, param, value):
    """
    Validate that the targets are valid (i.e. one of: all, assemble, annotate,
    qa or their respective sub-targets)
    """

    valid_targets = ['all',
                     'assemble', 'contig',
                     'annotate', 'annotate_genes', 'annotate_custom',
                     'qa', 'merqury', 'kraken2', 'fcs_gx', 'fcs_adaptor',
                     'mapping', 'pantools', 'busco', 'omark', 'mash', 'ntsynt',
                     'sans', 'pangrowth', 'statistics'
                     ]

    for target in value:
        if target not in valid_targets:
            raise click.BadParameter(f'Invalid target: {target}')

    return value


@cli.command('run',
    short_help='Run the MoGAAAP pipeline and create a report')
@click.option('--workdir', '-d',
    default='.',
    show_default=True,
    type=click.Path(exists=True),
    help='Working directory for MoGAAAP')
@click.option('--configfile', '-c',
    default='config/config.yaml',
    show_default=True,
    help='Configuration YAML file for MoGAAAP')
@click.option('--reportfile', '-r',
    default='report.html',
    show_default=True,
    help='Report file for MoGAAAP')
@click.option('--cores', '-j',
    default=multiprocessing.cpu_count(), # all available cores
    show_default=True,
    help='Number of cores to use')
@click.option('--memory', '-m',
    default=round(psutil.virtual_memory().total / (1024 ** 3) / 2), # half of total memory
    show_default=True,
    help='Amount of memory to use in GB (please put lower than what is available because this limit is not enforced)')
@click.option('--dryrun', '-n',
    is_flag=True,
    help='Dry run the pipeline')
@click.option('--other', '-o',
    help='Other options to pass to Snakemake')
@click.argument('targets',
    nargs=-1,
    required=False,
    callback=validate_targets)
def run(workdir, configfile, reportfile, cores, memory, dryrun, other, targets):
    """Run the MoGAAAP pipeline and create a report.
    You may select the target of the pipeline by specifying the target name(s) as arguments.
    Possible TARGETS are: all, assembly, scaffold, analyse, annotate, qa (default: all).
    """

    workdir = os.path.abspath(workdir)
    configfile = os.path.abspath(configfile)
    reportfile = os.path.abspath(reportfile)
    if not targets:
        targets = ["all"]

    run_mogaaap(workdir, configfile, reportfile, cores, memory, dryrun, other, targets)


def main():
    show_ascii_art()
    cli()


if __name__ == "__main__":
    main()
