rule get_reference:
    input:
        lambda wildcards: config["ref_genome"][wildcards.reference],
    output:
        "results/{asmname}/2.scaffolding/00.reference/{reference}.fa"
    log:
        "results/logs/2.scaffolding/get_reference/{asmname}/{reference}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/get_reference/{asmname}/{reference}.txt"
    shell:
        "ln -s {input} {output} &> {log}"
