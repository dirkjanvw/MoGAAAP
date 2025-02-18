rule kmerdb_build:
    input:
        lambda wildcards: expand("final_output/{asmname}.full.fa", asmname=get_all_accessions_from_asmset(wildcards.asmset)),
    output:
        "results/{asmset}/5.quality_assessment/08.kmer-db/{k}/{asmset}.db",
    log:
        "results/logs/5.quality_assessment/kmer-db_build/{k}/{asmset}.log"
    benchmark:
        "results/benchmarks/5.quality_assessment/kmer-db_build/{k}/{asmset}.txt"
    threads:
        10
    conda:
        "../../envs/kmer-db.yaml"
    shell:
        """(
        files=$(dirname {output})/files.list
        ls {input} > ${{files}}
        kmer-db build -k {wildcards.k} -t {threads} ${{files}} {output}
        ) &> {log}"""

rule kmerdb_compare:
    input:
        "results/{asmset}/5.quality_assessment/08.kmer-db/{k}/{asmset}.db",
    output:
        "results/{asmset}/5.quality_assessment/08.kmer-db/{k}/{asmset}.csv",
    log:
        "results/logs/5.quality_assessment/kmer-db_compare/{k}/{asmset}.log"
    benchmark:
        "results/benchmarks/5.quality_assessment/kmer-db_compare/{k}/{asmset}.txt"
    threads:
        10
    conda:
        "../../envs/kmer-db.yaml"
    shell:
        "kmer-db all2all {input} {output} &> {log}"

rule kmerdb_distance:
    input:
        "results/{asmset}/5.quality_assessment/08.kmer-db/{k}/{asmset}.csv",
    output:
        "results/{asmset}/5.quality_assessment/08.kmer-db/{k}/{asmset}.csv.jaccard",
        "results/{asmset}/5.quality_assessment/08.kmer-db/{k}/{asmset}.csv.min",
        "results/{asmset}/5.quality_assessment/08.kmer-db/{k}/{asmset}.csv.max",
        "results/{asmset}/5.quality_assessment/08.kmer-db/{k}/{asmset}.csv.cosine",
        "results/{asmset}/5.quality_assessment/08.kmer-db/{k}/{asmset}.csv.mash",
        "results/{asmset}/5.quality_assessment/08.kmer-db/{k}/{asmset}.csv.ani",
    log:
        "results/logs/5.quality_assessment/kmer-db_distance/{k}/{asmset}.log"
    benchmark:
        "results/benchmarks/5.quality_assessment/kmer-db_distance/{k}/{asmset}.txt"
    threads:
        10
    conda:
        "../../envs/kmer-db.yaml"
    shell:
        "kmer-db distance jaccard min max cosine mash ani {input} &> {log}"

rule kmerdb_heatmap:
    input:
        "results/{asmset}/5.quality_assessment/08.kmer-db/{k}/{asmset}.csv.mash",
    output:
        report("results/{asmset}/5.quality_assessment/08.kmer-db/{k}/{asmset}.csv.mash.pdf",
            category="Phylogeny",
            caption="../../report/kmer-db.rst",
            labels={"type": "kmer-db", "set": "{asmset}", "k": "{k}", "distance": "mash"}),
    log:
        "results/logs/5.quality_assessment/kmer-db_heatmap/{k}/{asmset}.log"
    benchmark:
        "results/benchmarks/5.quality_assessment/kmer-db_heatmap/{k}/{asmset}.txt"
    conda:
        "../../envs/rbase.yaml"
    shell:
        "Rscript workflow/scripts/pheatmap_lower_triangle.R {input} {output} &> {log}"
