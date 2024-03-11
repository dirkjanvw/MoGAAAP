# THIS SCRIPT SHOULD NOT BE RUN ON ITS OWN

def rename_sequences(scaffolds: str, unassigned: str, conversion_table: str,
                     chromosome_names: dict, output: str, prefix: str) -> None:
    """
    Rename the sequences in the scaffolds file using the conversion table.
    :param scaffolds: Filename of the scaffolds FASTA file
    :param unassigned: Filename of the unassigned FASTA file
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
    with open(scaffolds) as f1, open(unassigned) as f2, open(output, 'w') as out:
        for line in f1:
            if line.startswith('>'):
                scaffold = line.strip().lstrip('>')
                chrom = conversion.get(scaffold, scaffold)
                out.write(f'>{prefix}_Chr{chrom}\n')
            else:
                out.write(line)
        c = 1
        for line in f2:
            if line.startswith('>'):
                out.write(f'>{prefix}_Un{c}\n')
                c += 1
            else:
                out.write(line)


rename_sequences(str(snakemake.input.assigned), str(snakemake.input.unassigned),
        str(snakemake.input.table), snakemake.params.chroms, snakemake.output[0],
        snakemake.wildcards.asmname)
