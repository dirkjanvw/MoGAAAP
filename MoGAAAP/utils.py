import click


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
    click.secho('Initialising a new MoGAAAP pipeline', fg='blue')


def configure_mogaaap(workdir, configfile):
    click.secho('Configuring the MoGAAAP pipeline', fg='blue')


def run_mogaaap(workdir):
    click.secho('Running the MoGAAAP pipeline', fg='blue')
