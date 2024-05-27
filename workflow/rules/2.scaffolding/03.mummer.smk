rule get_reference:
    input:
        lambda wildcards: config["ref_genome"][wildcards.reference],
    output:
        temporary("results/{asmname}/2.scaffolding/03.mummer/{reference}.fa"),
    log:
        "results/logs/2.scaffolding/get_reference/{asmname}/{reference}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/get_reference/{asmname}/{reference}.txt"
    shell:
        "cp {input} {output} &> {log}"

rule mummer:
    input:
        reference = "results/{asmname}/2.scaffolding/03.mummer/{reference}.fa",
        assembly = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa",
    output:
        "results/{asmname}/2.scaffolding/03.mummer/{asmname}.vs.{reference}.delta",
    log:
        "results/logs/2.scaffolding/mummer/{reference}/{asmname}.vs.{reference}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/mummer/{reference}/{asmname}.vs.{reference}.txt"
    threads:
        24
    container:
        "workflow/singularity/mummer/mummer-4.0.0rc1.sif"
    shell:
        "nucmer -t {threads} -l 100 -g 100 --prefix=$(echo {output} | rev | cut -d '.' -f 2- | rev) {input.reference} {input.assembly} &> {log}"

rule dotplot:
    input:
        delta = "results/{asmname}/2.scaffolding/03.mummer/{asmname}.vs.{reference}.delta",
        reference = "results/{asmname}/2.scaffolding/03.mummer/{reference}.fa",
    output:
        gp = "results/{asmname}/2.scaffolding/03.mummer/{asmname}.vs.{reference}.plot.gp",
        png = report("results/{asmname}/2.scaffolding/03.mummer/{asmname}.vs.{reference}.plot.png",
            category="MUMmerplot",
            caption="../../report/mummerplot.rst",
            labels={"assembly": "{asmname}",
                    "stage": "scaffolds"}),
    log:
        "results/logs/2.scaffolding/dotplot/{asmname}/{asmname}.vs.{reference}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/dotplot/{asmname}/{asmname}.vs.{reference}.txt"
    container:
        "workflow/singularity/mummer/mummer-4.0.0rc1.sif"
    shell:  #assumes input and output are in same directory
        """
        (
        cd $(dirname {output.gp})
        mummerplot --filter --png --large --prefix=$(basename {output.gp} | rev | cut -d'.' -f 2- | rev) --title {wildcards.asmname} $(basename {input.delta})
        ) &> {log}
        """
