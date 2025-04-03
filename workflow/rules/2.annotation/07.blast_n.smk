def get_query_file_blast_n(wildcards):
    if config["organellar"].get(wildcards.query_name):
        query_file = config["organellar"][wildcards.query_name]
    else:
        query_file = config["nucl_queries"][wildcards.query_name]

    return query_file

rule blast_n:
    input:
        query_file = get_query_file_blast_n,
        blast_db = "results/{asmname}/2.annotation/05.blastdb/{asmname}.BDB",
    output:
        "results/{asmname}/2.annotation/07.blast_n/{query_name}.vs.{asmname}.m7"
    log:
        "results/logs/2.annotation/blast_n/{asmname}/{query_name}.vs.{asmname}.log"
    benchmark:
        "results/benchmarks/2.annotation/blast_n/{asmname}/{query_name}.vs.{asmname}.txt"
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

rule blast_n_to_tsv:
    input:
        "results/{asmname}/2.annotation/07.blast_n/{query_name}.vs.{asmname}.m7",
    output:
        "results/{asmname}/2.annotation/07.blast_n/{query_name}.vs.{asmname}.tsv",
    log:
        "results/logs/2.annotation/blast_n_to_tsv/{asmname}/{query_name}.vs.{asmname}.log"
    benchmark:
        "results/benchmarks/2.annotation/blast_n_to_tsv/{asmname}/{query_name}.vs.{asmname}.txt"
    shell:
        "awk 'BEGIN{{FS = OFS = \"\\t\";}} /^# Fields: /{{$1=substr($1,11); gsub(/, /,\"\\t\",$1); print $1; next;}} /^#/{{next;}} {{print;}}' {input} > {output} 2> {log}"

rule visualise_blast_n:
    input:
        "results/{asmname}/2.annotation/07.blast_n/{query_name}.vs.{asmname}.tsv",
    output:
        report("results/{asmname}/2.annotation/07.blast_n/{query_name}.vs.{asmname}.html",
            category="Analysis",
            caption="../../report/blast.rst",
            labels={"asmname": "{asmname}", "query_name": "{query_name}"}),
    log:
        "results/logs/2.annotation/visualise_blast_n/{asmname}/{query_name}.vs.{asmname}.log"
    benchmark:
        "results/benchmarks/2.annotation/visualise_blast_n/{asmname}/{query_name}.vs.{asmname}.txt"
    conda:
        "../../envs/csvtotable.yaml"
    shell:
        "csvtotable -d $'\\t' {input} {output} &> {log}"
