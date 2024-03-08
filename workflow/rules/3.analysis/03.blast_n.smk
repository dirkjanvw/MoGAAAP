rule blast_n:
    input:
        query_file = lambda wildcards: config["nucl_queries"][wildcards.query_name],
        blast_db = "results/{asmname}/3.analysis/01.blastdb/{asmname}.BDB",
    output:
        "results/{asmname}/3.analysis/03.blast_n/{query_name}.vs.{asmname}.m7"
    log:
        "results/logs/blast_n/{asmname}/{query_name}.vs.{asmname}.log"
    benchmark:
        "results/benchmarks/3.analysis/blast_n/{asmname}/{query_name}.vs.{asmname}.txt"
    threads:
        12
    conda:
        "../../envs/blast.yaml"
    shell:
        """
        blastn -task blastn \
           -query {input.query_file} -db {input.blast_db} -out {output} \
           -outfmt '7 qseqid sseqid evalue pident score length nident mismatch gaps frames qstart qend sstart send qcovhsp qlen slen' \
           -evalue 1e-120 -dbsize 1000000 -dust no -word_size 24 -xdrop_ungap 50 -xdrop_gap 500 -xdrop_gap_final 1000 \
           -max_target_seqs 12000 -num_threads {threads} &> {log}
        """
