def get_blast_db(wildcards): 
    return ["results/01.blastdb/{asmname}/{asmname}.BDB"]

rule blast_n:
    input:
        query_file = lambda wildcards: config["nucl_queries"][wildcards.query_name], 
        blast_db   = get_blast_db, 
    output:
        "results/03.blast_n/{asmname}/{query_name}.vs.{asmname}.m7" 
    log:
        "results/logs/03.blast_n/{asmname}/{query_name}.vs.{asmname}.log" 
    benchmark:
        "results/benchmarks/03.blast_n/{asmname}/{query_name}.vs.{asmname}.txt" 
    threads:
        12 
    shell:
        """
        /path_to_data/x-bin-x/ncbi-blast-2.15/bin/blastn -task blastn \
           -query {input.query_file} -db {input.blast_db} -out {output} \
           -outfmt '7 qseqid sseqid evalue pident score length nident mismatch gaps frames qstart qend sstart send qcovhsp qlen slen' \
           -evalue 1e-120 -dbsize 1000000 -dust no -word_size 24 -xdrop_ungap 50 -xdrop_gap 500 -xdrop_gap_final 1000 \
           -max_target_seqs 12000 -num_threads {threads} &> {log} 
        """
