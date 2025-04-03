rule blast_p:
    input:
        query_file = lambda wildcards: config["prot_queries"][wildcards.query_name],
        blast_db   = "results/{asmname}/3.analysis/01.blastdb/{asmname}.BDB",
    output:
        "results/{asmname}/3.analysis/02.blast_p/{query_name}.vs.{asmname}.m7"
    log:
        "results/logs/3.analysis/blast_p/{asmname}/{query_name}.vs.{asmname}.log"
    benchmark:
        "results/benchmarks/3.analysis/blast_p/{asmname}/{query_name}.vs.{asmname}.txt"
    threads:
        12
    conda:
        "../../envs/blast.yaml"
    shell:
        """
        tblastn \
           -query {input.query_file} -db {input.blast_db} -out {output} \
           -outfmt '7 qseqid sseqid evalue pident score length nident mismatch gaps frames qstart qend sstart send qcovhsp qlen slen' \
           -evalue 1e-10 -dbsize 1000000 -seg no -word_size 3 -max_intron_length 120 -qcov_hsp_perc 20 -xdrop_ungap 40 -xdrop_gap 60 \
           -xdrop_gap_final 200 -max_target_seqs 12000 -num_threads {threads} &> {log}
        """

rule blast_p_to_tsv:
    input:
        "results/{asmname}/3.analysis/02.blast_p/{query_name}.vs.{asmname}.m7",
    output:
        "results/{asmname}/3.analysis/02.blast_p/{query_name}.vs.{asmname}.tsv",
    log:
        "results/logs/3.analysis/blast_p_to_tsv/{asmname}/{query_name}.vs.{asmname}.log"
    benchmark:
        "results/benchmarks/3.analysis/blast_p_to_tsv/{asmname}/{query_name}.vs.{asmname}.txt"
    shell:
        "awk 'BEGIN{{FS = OFS = \"\\t\";}} /^# Fields: /{{$1=substr($1,11); gsub(/, /,\"\\t\",$1); print $1; next;}} /^#/{{next;}} {{print;}}' {input} > {output} 2> {log}"

rule visualise_blast_p:
    input:
        "results/{asmname}/3.analysis/02.blast_p/{query_name}.vs.{asmname}.tsv",
    output:
        report("results/{asmname}/3.analysis/02.blast_p/{query_name}.vs.{asmname}.html",
            category="Analysis",
            caption="../../report/blast.rst",
            labels={"asmname": "{asmname}", "query_name": "{query_name}"}),
    log:
        "results/logs/3.analysis/visualise_blast_p/{asmname}/{query_name}.vs.{asmname}.log"
    benchmark:
        "results/benchmarks/3.analysis/visualise_blast_p/{asmname}/{query_name}.vs.{asmname}.txt"
    conda:
        "../../envs/csvtotable.yaml"
    shell:
        "csvtotable -d $'\\t' {input} {output} &> {log}"
