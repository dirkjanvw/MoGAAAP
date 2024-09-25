# THIS SCRIPT SHOULD NOT BE RUN ON ITS OWN
import math

def rename_sequences(sequences: str, conversion_table: str,
                     chromosome_names: dict, output: str, prefix: str) -> None:
    """
    Rename the sequences in the input file using the conversion table.
    :param sequences: Filename of the sequences to be renamed
    :param conversion_table: Filename of the conversion table (first column:
        reference name, second column: scaffolds name)
    :param chromosome_names: Dictionary of new names for the chromosomes (first
        column: chromosome number, second column: reference name)
    :param output: Filename of the output FASTA file
    :param prefix: Prefix for the chromosome names
    :return: None
    """

    # Read the conversion table
    conversion = {}
    with open(conversion_table) as f:
        for line in f:
            ref, scaffold = line.strip().split()
            for k, v in chromosome_names.items():
                if ref == v:
                    conversion[scaffold] = k
                    break

    # Read the input file and write the output
    sequencedict = {}
    with open(sequences) as infile:
        for line in infile:
            if line.startswith('>'):
                sequencename = line.strip().lstrip('>')
                sequencedict[sequencename] = ""
            else:
                sequencedict[sequencename] += line

    # Write the output
    chromosomeleadingzeroes = int(math.log10(len(conversion)))+1
    contigleadingzeroes = int(math.log10(len(sequencedict)))+1
    c = 1  #for counting unnamed sequences
    with open(output, 'w') as outfile:
        for sequencename, sequence in sequencedict.items():
            if sequencename in conversion:
                outfile.write(f'>{prefix}_Chr{str(conversion.get(sequencename, sequencename)).zfill(chromosomeleadingzeroes)}\n')
            else:
                outfile.write(f'>{prefix}_Un{str(c).zfill(contigleadingzeroes)} {sequencename}\n')
                c += 1
            outfile.write(sequence)


rename_sequences(str(snakemake.input.all),
                 str(snakemake.input.table), snakemake.params.chroms,
                 snakemake.output[0],
                 snakemake.wildcards.asmname)
