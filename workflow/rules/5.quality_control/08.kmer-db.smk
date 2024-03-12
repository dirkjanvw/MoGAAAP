rule kmerdb_build:
    input:
        lambda wildcards: expand("results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa", asmname=config["set"][wildcards.asmset]),
    output:
        "results/{asmset}/5.quality_control/08.kmer-db/k{k}.db",
    log:
        "results/logs/5.quality_control/kmer-db_build/{k}/{asmset}.log"
    benchmark:
        "results/benchmarks/5.quality_control/kmer-db_build/{k}/{asmset}.txt"
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
        "results/{asmset}/5.quality_control/08.kmer-db/k{k}.db",
    output:
        "results/{asmset}/5.quality_control/08.kmer-db/k{k}.csv",
    log:
        "results/logs/5.quality_control/kmer-db_compare/{k}/{asmset}.log"
    benchmark:
        "results/benchmarks/5.quality_control/kmer-db_compare/{k}/{asmset}.txt"
    threads:
        10
    conda:
        "../../envs/kmer-db.yaml"
    shell:
        "kmer-db all2all {input} {output} &> {log}"

rule kmerdb_distance:
    input:
        "results/{asmset}/5.quality_control/08.kmer-db/k{k}.csv",
    output:
        "results/{asmset}/5.quality_control/08.kmer-db/k{k}.csv.jaccard",
        "results/{asmset}/5.quality_control/08.kmer-db/k{k}.csv.min",
        "results/{asmset}/5.quality_control/08.kmer-db/k{k}.csv.max",
        "results/{asmset}/5.quality_control/08.kmer-db/k{k}.csv.cosine",
        "results/{asmset}/5.quality_control/08.kmer-db/k{k}.csv.mash",
        "results/{asmset}/5.quality_control/08.kmer-db/k{k}.csv.ani",
    log:
        "results/logs/5.quality_control/kmer-db_distance/{k}/{asmset}.log"
    benchmark:
        "results/benchmarks/5.quality_control/kmer-db_distance/{k}/{asmset}.txt"
    threads:
        10
    conda:
        "../../envs/kmer-db.yaml"
    shell:
        "kmer-db distance jaccard min max cosine mash ani {input} &> {log}"

rule kmerdb_heatmap:
    input:
        "results/{asmset}/5.quality_control/08.kmer-db/k{k}.csv.mash",
    output:
        report("results/{asmset}/5.quality_control/08.kmer-db/k{k}.csv.mash.pdf", category="General", labels={"type": "kmer-db", "set": "{asmset}", "k": "{k}", "distance": "mash"}),
    log:
        "results/logs/5.quality_control/kmer-db_heatmap/{k}/{asmset}.log"
    benchmark:
        "results/benchmarks/5.quality_control/kmer-db_heatmap/{k}/{asmset}.txt"
    conda:
        "../../envs/rbase.yaml"
    shell:
        "Rscript workflow/scripts/pheatmap_lower_triangle.R {input} {output} &> {log}"