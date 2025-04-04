rule blastdb:
    input:
        "final_output/{asmname}.full.fa",
    output:
        "results/{asmname}/2.annotation/05.blastdb/{asmname}.BDB",
    log:
        "results/logs/2.annotation/blastdb/{asmname}.log"
    benchmark:
        "results/benchmarks/2.annotation/blastdb/{asmname}.txt"
    conda:
        "../../envs/blast.yaml"
    shell:
        """
        makeblastdb -in {input} -out {output} \
          -dbtype nucl -input_type fasta -max_file_sz 2GB -hash_index -parse_seqids &> {log}
        touch {output}
        date >> {output}
        """
