rule mummer:
    input:
        reference = lambda wildcards: config["ref_genome"][wildcards.reference],
        assembly = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa",
    output:
        "results/{asmname}/2.scaffolding/03.mummer/{asmname}.vs.{reference}.delta",
    log:
        "results/logs/2.scaffolding/mummer/{reference}/{asmname}.vs.{reference}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/mummer/{reference}/{asmname}.vs.{reference}.txt"
    threads:
        24
    conda:
        "../../envs/mummer.yaml"
    shell:
        "nucmer -t {threads} -l 1000 -g 1000 --prefix=$(echo {output} | rev | cut -d '.' -f 2- | rev) {input.reference} {input.assembly} &> {log}"

rule dotplot:
    input:
        "results/{asmname}/2.scaffolding/03.mummer/{asmname}.vs.{reference}.delta",
    output:
        "results/{asmname}/2.scaffolding/03.mummer/{asmname}.vs.{reference}.plot.gp",
    log:
        "results/logs/2.scaffolding/dotplot/{asmname}/{asmname}.vs.{reference}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/dotplot/{asmname}/{asmname}.vs.{reference}.txt"
    conda:
        "../../envs/mummer.yaml"
    shell:  #assumes input and output are in same directory
        """
        (
        cd $(dirname {output})
        mummerplot --filter --png --large --prefix=$(basename {output} | rev | cut -d'.' -f 2- | rev) --title {wildcards.asmname} $(basename {input})
        ) &> {log}
        """
