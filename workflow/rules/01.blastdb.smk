rule blastdb:
    input:
        lambda wildcards: config["assemblies"][wildcards.asmname]
    output:
        "results/01.blastdb/{asmname}/{asmname}.BDB"
    log:
        "results/logs/01.blastdb/{asmname}.log"
    benchmark:
        "results/benchmarks/01.blastdb/{asmname}.txt"
    shell:
        """
        /path_to_data/x-bin-x/ncbi-blast-2.15/bin/makeblastdb -in {input} -out {output} \
          -dbtype nucl -input_type fasta -max_file_sz 2GB -hash_index -parse_seqids &> {log} 
        touch {output} 
        date >> {output}
        """
