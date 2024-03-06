def get_blast_db(wildcards): 
    return ["results/01.blastdb/{asmname}/{asmname}.BDB"]

rule blast_p:
    input:
        query_file = lambda wildcards: config["prot_queries"][wildcards.query_name], 
        blast_db   = get_blast_db, 
    output:
        "results/02.blast_p/{asmname}/{query_name}.vs.{asmname}.m7" 
    log:
        "results/logs/02.blast_p/{asmname}/{query_name}.vs.{asmname}.log" 
    benchmark:
        "results/benchmarks/02.blast_p/{asmname}/{query_name}.vs.{asmname}.txt" 
    threads:
        12 
    shell:
        """
        /path_to_data/x-bin-x/ncbi-blast-2.15/bin/tblastn \
           -query {input.query_file} -db {input.blast_db} -out {output} \
           -outfmt '7 qseqid sseqid evalue pident score length nident mismatch gaps frames qstart qend sstart send qcovhsp qlen slen' \
           -evalue 1e-10 -dbsize 1000000 -seg no -word_size 3 -max_intron_length 120 -qcov_hsp_perc 20 -xdrop_ungap 40 -xdrop_gap 60 \
           -xdrop_gap_final 200 -max_target_seqs 12000 -num_threads {threads} &> {log} 
        """
