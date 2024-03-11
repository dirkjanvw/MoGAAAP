rule create_agat_config:
    output:
        "results/{asmname}/agat_config.yaml"
    log:
        "results/logs/create_agat_config/{asmname}.log"
    benchmark:
        "results/benchmarks/create_agat_config/{asmname}.txt"
    conda:
        "../envs/agat.yaml"
    shell:
        "agat config --expose --output {output} --no-log --no-progress_bar --tabix --prefix_new_id nbis_{wildcards.asmname} &> {log}"
