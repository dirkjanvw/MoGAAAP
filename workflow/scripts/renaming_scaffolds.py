# THIS SCRIPT SHOULD NOT BE RUN ON ITS OWN

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
    c = 1  #for counting unnamed sequences
    with open(sequences) as f,\
            open(output, 'w') as out:
        for line in f:
            if line.startswith('>'):
                sequence = line.strip().lstrip('>')
                if sequence in conversion:
                    out.write(f'>{prefix}_Chr{conversion.get(sequence, sequence)}\n')  #TODO: Add leading zeroes
                else:
                    out.write(f'>{prefix}_Un{c}\n')  #TODO: Add leading zeroes
                    c += 1
            else:
                out.write(line)


rename_sequences(str(snakemake.input.all),
                 str(snakemake.input.table), snakemake.params.chroms,
                 snakemake.output[0],
                 snakemake.wildcards.asmname)
