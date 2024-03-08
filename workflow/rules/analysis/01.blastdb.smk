rule blastdb:
    input:
        "results/{asmname}/2.scaffolding/03.renaming/{asmname}.fa",
    output:
        "results/{asmname}/3.analysis/01.blastdb/{asmname}.BDB"
    log:
        "results/logs/01.blastdb/{asmname}.log"
    benchmark:
        "results/benchmarks/01.blastdb/{asmname}.txt"
    conda:
        "../../envs/blast.yaml"
    shell:
        """
        makeblastdb -in {input} -out {output} \
          -dbtype nucl -input_type fasta -max_file_sz 2GB -hash_index -parse_seqids &> {log}
        touch {output}
        date >> {output}
        """
