#!/usr/bin/env python3
# Written by Dirk-Jan van Workum (2025)
# This is a wrapper around the MoGAAAP Snakemake workflow
#
# Usage: mogaaap [subcommand] [options]
#
# Subcommands:
# - init: Initialise a new MoGAAAP pipeline
# - configure: Configure the MoGAAAP pipeline
# - run: Run the MoGAAAP pipeline and create a report

import click
from importlib.metadata import version
import multiprocessing, psutil, os
from .utils import show_ascii_art, init_mogaaap, configure_mogaaap, run_mogaaap


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


@click.group(cls=OrderedGroup, context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(version=version("MoGAAAP"))
def cli():
    """This is a wrapper script around the MoGAAAP Snakemake workflow."""
    pass


@cli.command('init',
    short_help='Initialise a new MoGAAAP pipeline')
@click.option('--workdir', '-d',
    default='.',
    show_default=True,
    help='Working directory for MoGAAAP')
def init(workdir):
    """Initialise a new MoGAAAP pipeline"""

    workdir = os.path.abspath(workdir)

    init_mogaaap(workdir)


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
    required=True,
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
    required=True,
    help='ODB name for BUSCO')
@click.option('--helixer-model', '-e',
    required=True,
    help='Helixer model name (land_plant, vertebrate, invertebrate, fungi)')
@click.option('--gxdb', '-x',
    required=True,
    type=click.Path(exists=True),
    help='GXDB database location')
@click.option('--omadb', '-o',
    required=True,
    type=click.Path(exists=True),
    help='OMA database location')
@click.option('--kraken2db', '-k',
    required=True,
    type=click.Path(exists=True),
    help='kraken2 database location')
def configure(workdir, samples, reference_fasta, reference_gff, mitochondrion,
    chloroplast, telomere, odb, helixer_model, gxdb, omadb, kraken2db):
    """Configure the MoGAAAP pipeline"""

    workdir = os.path.abspath(workdir)
    samples = os.path.abspath(samples)
    reference_fasta = os.path.abspath(reference_fasta)
    reference_gff = os.path.abspath(reference_gff)
    mitochondrion = os.path.abspath(mitochondrion)
    if chloroplast:
        chloroplast = os.path.abspath(chloroplast)
    gxdb = os.path.abspath(gxdb)
    omadb = os.path.abspath(omadb)
    kraken2db = os.path.abspath(kraken2db)

    configure_mogaaap(workdir, samples, reference_fasta, reference_gff,
        mitochondrion, chloroplast, telomere, odb, helixer_model, gxdb, omadb,
        kraken2db)


def validate_targets(ctx, param, value):
    """
    Validate that the targets are valid (i.e. one of: all, assembly, scaffold,
        analyse, annotate, qa)
    """

    valid_targets = ['all', 'assemble', 'contig', 'annotate', 'annotate_genes',
        'annotate_custom', 'qa']

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
    help='Amount of memory to use in GB')
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

    run_mogaaap(workdir, configfile, reportfile, cores, memory, dryrun, other, targets)


def main():
    show_ascii_art()
    cli()


if __name__ == "__main__":
    main()
