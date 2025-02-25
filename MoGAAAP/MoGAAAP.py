#!/usr/bin/env python3
# Written by Dirk-Jan van Workum (2025)
# This is a wrapper around the MoGAAAP Snakemake workflow
#
# Usage: mogaaap [subcommand] [options]
#
# Subcommands:
# - init: Initialise a new MoGAAAP pipeline
# - config: Configure the MoGAAAP pipeline
# - run: Run the MoGAAAP pipeline and create a report

import click
from importlib.metadata import version


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(version=version("MoGAAAP"))
def cli():
    """This is a wrapper script around the MoGAAAP Snakemake workflow"""
    pass


@cli.command('init',
    short_help='Initialise a new MoGAAAP pipeline')
@click.option('--workdir', '-d',
    default='.',
    help='Working directory for MoGAAAP')
def init(workdir):
    """Initialise a new MoGAAAP pipeline"""

    click.secho('Initialising a new MoGAAAP pipeline', fg='blue')

    #TODO: Implement init function


@cli.command('config',
    short_help='Configure the MoGAAAP pipeline')
@click.option('--workdir', '-d',
    default='.',
    help='Working directory for MoGAAAP')
@click.option('--configfile', '-c',
    default='config/config.yaml',
    help='Configuration file for MoGAAAP')
def config(workdir, configfile):
    """Configure the MoGAAAP pipeline"""

    click.secho('Configuring the MoGAAAP pipeline', fg='blue')

    #TODO: Implement config function


@cli.command('run',
    short_help='Run the MoGAAAP pipeline and create a report')
@click.option('--workdir', '-d',
    default='.',
    help='Working directory for MoGAAAP')
def run(workdir):
    """Run the MoGAAAP pipeline and create a report"""

    click.secho('Running the MoGAAAP pipeline', fg='blue')


def main():
    #show_ascii_art() #TODO
    cli()


if __name__ == "__main__":
    main()
