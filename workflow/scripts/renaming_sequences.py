def rename_sequences(scaffolds: str, conversion_table: str,
                     chromosome_names: dict, output: str, prefix: str) -> None:
    """
    Rename the sequences in the scaffolds file using the conversion table.
    :param scaffolds: Filename of the scaffolds FASTA file
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

    # Read the scaffolds file and write the output
    with open(scaffolds) as f, open(output, 'w') as out:
        for line in f:
            if line.startswith('>'):
                scaffold = line.strip().lstrip('>')
                chrom = conversion.get(scaffold, scaffold)
                out.write(f'>{prefix}_Chr{chrom}\n')
            else:
                out.write(line)


rename_sequences(snakemake.input[0], snakemake.input[1], snakemake.params[0],
                 snakemake.output[0], snakemake.wildcards.asmname)
