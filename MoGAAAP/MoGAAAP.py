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
    help='Working directory for MoGAAAP')
def init(workdir):
    """Initialise a new MoGAAAP pipeline"""

    init_mogaaap(workdir)


@cli.command('configure',
    short_help='Configure the MoGAAAP pipeline')
@click.option('--workdir', '-d',
    default='.',
    help='Working directory for MoGAAAP')
@click.option('--configfile', '-c',
    default='config/config.yaml',
    help='Configuration file for MoGAAAP')
def configure(workdir, configfile):
    """Configure the MoGAAAP pipeline"""

    configure_mogaaap(workdir, configfile)


@cli.command('run',
    short_help='Run the MoGAAAP pipeline and create a report')
@click.option('--workdir', '-d',
    default='.',
    help='Working directory for MoGAAAP')
def run(workdir):
    """Run the MoGAAAP pipeline and create a report"""

    run_mogaaap(workdir)


def main():
    show_ascii_art()
    cli()


if __name__ == "__main__":
    main()
