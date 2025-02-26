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
import multiprocessing, psutil
from .utils import show_ascii_art, init_mogaaap, configure_mogaaap, run_mogaaap


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
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

    init_mogaaap(workdir)


@cli.command('configure',
    short_help='Configure the MoGAAAP pipeline')
@click.option('--workdir', '-d',
    default='.',
    show_default=True,
    help='Working directory for MoGAAAP')
@click.option('--configfile', '-c',
    default='config/config.yaml',
    show_default=True,
    help='Configuration file for MoGAAAP')
def configure(workdir, configfile):
    """Configure the MoGAAAP pipeline"""

    configure_mogaaap(workdir, configfile)


def validate_targets(ctx, param, value):
    """
    Validate that the targets are valid (i.e. one of: all, assembly, scaffold,
        analyse, annotate, qa)
    """

    valid_targets = ['all', 'assembly', 'scaffold', 'analyse', 'annotate', 'qa']

    for target in value:
        if target not in valid_targets:
            raise click.BadParameter(f'Invalid target: {target}')

    return value


@cli.command('run')
@click.option('--workdir', '-d',
    default='.',
    show_default=True,
    help='Working directory for MoGAAAP')
@click.option('--configfile', '-c',
    default='config/config.yaml',
    show_default=True,
    type=click.Path(exists=True),
    help='Configuration YAML file for MoGAAAP')
@click.option('--reportfile', '-r',
    default='report.html',
    show_default=True,
    help='Report file for MoGAAAP')
@click.option('--cores', '-j',
    default=multiprocessing.cpu_count(),
    show_default=True,
    help='Number of cores to use')
@click.option('--memory', '-m',
    default=round(psutil.virtual_memory().total / (1024 ** 3)),
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

    run_mogaaap(workdir, configfile, reportfile, cores, memory, dryrun, other, targets)


def main():
    show_ascii_art()
    cli()


if __name__ == "__main__":
    main()
