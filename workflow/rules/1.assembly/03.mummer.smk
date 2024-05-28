rule get_reference_contigs:
    input:
        lambda wildcards: config["reference_genomes"][wildcards.reference]["genome"]
    output:
        temporary("results/{asmname}/1.assembly/03.mummer/{reference}.fa"),
    log:
        "results/logs/1.assembly/get_reference/{asmname}/{reference}.log"
    benchmark:
        "results/benchmarks/1.assembly/get_reference/{asmname}/{reference}.txt"
    shell:
        "cp {input} {output} &> {log}"

rule mummer_contigs:
    input:
        reference = "results/{asmname}/1.assembly/03.mummer/{reference}.fa",
        assembly = "results/{asmname}/1.assembly/02.contigs/{asmname}.min{minlen}.sorted.renamed.fa",
    output:
        "results/{asmname}/1.assembly/03.mummer/{asmname}.min{minlen}.vs.{reference}.delta",
    log:
        "results/logs/1.assembly/mummer/{reference}/{asmname}.min{minlen}.vs.{reference}.log"
    benchmark:
        "results/benchmarks/1.assembly/mummer/{reference}/{asmname}.min{minlen}.vs.{reference}.txt"
    threads:
        24
    container:
        "workflow/singularity/mummer/mummer-4.0.0rc1.sif"
    shell:
        "nucmer -t {threads} -l 1000 -g 1000 --prefix=$(echo {output} | rev | cut -d '.' -f 2- | rev) {input.reference} {input.assembly} &> {log}"

rule dotplot_contigs:
    input:
        delta = "results/{asmname}/1.assembly/03.mummer/{asmname}.min{minlen}.vs.{reference}.delta",
        reference = "results/{asmname}/1.assembly/03.mummer/{reference}.fa",
    output:
        gp = "results/{asmname}/1.assembly/03.mummer/{asmname}.min{minlen}.vs.{reference}.plot.gp",
        png = report("results/{asmname}/1.assembly/03.mummer/{asmname}.min{minlen}.vs.{reference}.plot.png",
            category="MUMmerplot",
            caption="../../report/mummerplot.rst",
            labels={"assembly": "{asmname}",
                    "stage": "contigs"}),
    log:
        "results/logs/1.assembly/dotplot/{asmname}/{asmname}.min{minlen}.vs.{reference}.log"
    benchmark:
        "results/benchmarks/1.assembly/dotplot/{asmname}/{asmname}.min{minlen}.vs.{reference}.txt"
    container:
        "workflow/singularity/mummer/mummer-4.0.0rc1.sif"
    shell:  #assumes input and output are in same directory
        """
        (
        cd $(dirname {output.gp})
        mummerplot --filter --png --large --prefix=$(basename {output.gp} | rev | cut -d'.' -f 2- | rev) --title {wildcards.asmname} $(basename {input.delta})
        ) &> {log}
        """
