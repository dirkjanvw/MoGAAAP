#!/usr/bin/env python3

import argparse
import math


def parse_args() -> argparse.Namespace:
    """
    Parse input arguments (karyotype file, other input files and output file)
    :return: argparse.Namespace
    """
    parser = argparse.ArgumentParser(description="Create a circos "
                                                 "configuration file")
    parser.add_argument("-k", "--karyotype", type=str,
                        help="Karyotype file",
                        required=True)
    parser.add_argument("-o", "--output", type=str,
                        help="Configuration file (default: circos.conf)",
                        default="circos.conf")
    parser.add_argument("circos_files", nargs="+",
                        help="Other input files for circos; max for .fract.circos files in 1, for all other files it is 100")
    return parser.parse_args()


def parse_sequence_names(karyotype: str) -> list:
    """
    Parse sequence names from a karyotype file
    :param karyotype: str
    :return: list
    """
    sequences = []
    with open(karyotype, "r") as k:
        for line in k:
            if line.startswith("chr"):
                sequences.append(line.split()[3])
    return sequences


def write_circos_config(karyotype: str, circos_files: list, output: str,
                        lower_bound: float = 0.5,
                        upper_bound: float = 0.99) -> None:
    """
    Write a circos configuration file
    :param karyotype: str
    :param circos_files: list
    :param output: str
    :param lower_bound: float
    :param upper_bound: float
    :return: None
    """

    # Parse sequence names from karyotype file
    sequences = parse_sequence_names(karyotype)

    # Add header
    header = f"""# circos.conf

# LIST OF CHROMOSOMES AND THEIR LENGTH #
karyotype = {karyotype}
chromosomes_units   = 1000000
chromosomes_display_default = yes
# CHROMOSOME RULE LABELS #
<<include ticks.conf>>

<ideogram>
# CHROMOSOME DISPLAY #
<spacing>
default = 10u
<pairwise {sequences[0]},{sequences[-1]}>
spacing = 60u
</pairwise>
</spacing>
radius    = 0.8r
thickness = 30p
fill      = yes
show_label = yes
label_font = bold
label_radius = 1r + 90p
label_size = 24
label_parallel = yes
# END OF CHROMOSOMES #
</ideogram>

<image>
angle_offset* = -87
</image>
"""

    # Add plots for each circos file
    steps = (upper_bound - lower_bound) / len(circos_files)
    plots = """
### PLOTS SECTION ###
<plots>
"""
    for i, circos_file in enumerate(circos_files):
        low = round(lower_bound + i * steps, 2)
        high = math.ceil((lower_bound + (i + 1) * steps) * 100) / 100
        maximum = 1 if circos_file.endswith(".fract.circos") else 100
        plots += f"""
# {circos_file} #
<plot>
type        = histogram
file        = {circos_file}
r1          = {high}r
r0          = {low}r
orientation = out
color       = black
fill_color  = lgrey
extend_bin  = no
min         = 0
max         = {maximum}
<backgrounds>
<background>
color       = vvlgrey
</background>
</backgrounds>
</plot>
"""
    plots += """</plots>\n"""

    # Add footer
    footer = f"""
################################################################
# The remaining content is standard and required. It is imported
# from default files in the Circos distribution.
#
# These should be present in every Circos configuration file and
# overridden as required. To see the content of these files,
# look in etc/ in the Circos distribution.

<image>
# Included from Circos distribution.
<<include etc/image.conf>>
</image>

# RGB/HSV color definitions, color lists, location of fonts, fill patterns.
# Included from Circos distribution.
<<include etc/colors_fonts_patterns.conf>>

# Debugging, I/O and other system parameters
# Included from Circos distribution.
<<include etc/housekeeping.conf>>
################################################################
"""

    with open(output, "w") as out:
        out.write(header)
        out.write(plots)
        out.write(footer)


def main():
    # Parse input arguments
    args = parse_args()

    # Write circos configuration file
    write_circos_config(args.karyotype, args.circos_files, args.output)


if __name__ == "__main__":
    main()
