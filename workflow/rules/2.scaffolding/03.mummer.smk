rule get_reference:
    input:
        lambda wildcards: config["reference_genomes"][wildcards.reference]["genome"],
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
    params:
        nucmer_maxgap = config["nucmer_maxgap"],
        nucmer_minmatch = config["nucmer_minmatch"],
    threads:
        lambda wildcards: len(config["reference_genomes"][get_reference_id(wildcards.asmname)]["chromosomes"]),  #the number of chromosomes
    container:
        "workflow/singularity/mummer/mummer-4.0.0rc1.sif"
    shell:
        "nucmer -t {threads} -l {params.nucmer_maxgap} -g {params.nucmer_minmatch} --prefix=$(echo {output} | rev | cut -d '.' -f 2- | rev) {input.reference} {input.assembly} &> {log}"

rule dotplot:
    input:
        delta = "results/{asmname}/2.scaffolding/03.mummer/{asmname}.vs.{reference}.delta",
        reference = "results/{asmname}/2.scaffolding/03.mummer/{reference}.fa",
    output:
        gp = "results/{asmname}/2.scaffolding/03.mummer/{asmname}.vs.{reference}.plot.gp",
        rplot = "results/{asmname}/2.scaffolding/03.mummer/{asmname}.vs.{reference}.plot.rplot",
        fplot = "results/{asmname}/2.scaffolding/03.mummer/{asmname}.vs.{reference}.plot.fplot",
        filterfile = "results/{asmname}/2.scaffolding/03.mummer/{asmname}.vs.{reference}.plot.filter",
        png = report("results/{asmname}/2.scaffolding/03.mummer/{asmname}.vs.{reference}.plot.png",
            category="MUMmerplot",
            caption="../../report/mummerplot.rst",
            labels={"assembly": "{asmname}",
                    "size": "default",
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

rule dotplot_large:
    input:
        gp = "results/{asmname}/2.scaffolding/03.mummer/{asmname}.vs.{reference}.plot.gp",
        rplot = "results/{asmname}/2.scaffolding/03.mummer/{asmname}.vs.{reference}.plot.rplot",
        fplot = "results/{asmname}/2.scaffolding/03.mummer/{asmname}.vs.{reference}.plot.fplot",
    output:
        gp = "results/{asmname}/2.scaffolding/03.mummer/{asmname}.vs.{reference}.plot.large.gp",
        png = report("results/{asmname}/2.scaffolding/03.mummer/{asmname}.vs.{reference}.plot.large.png",
            category="MUMmerplot",
            caption="../../report/mummerplot.rst",
            labels={"assembly": "{asmname}",
                    "size": "large",
                    "stage": "scaffolds"}),
    log:
        "results/logs/2.scaffolding/dotplot/{asmname}/{asmname}.vs.{reference}.large.log"
    benchmark:
        "results/benchmarks/2.scaffolding/dotplot/{asmname}/{asmname}.vs.{reference}.large.txt"
    container:
        "workflow/singularity/mummer/mummer-4.0.0rc1.sif"
    shell:
        """
        (
        sed 's/set terminal png tiny size 1400,1400/set terminal png large size 2600,2600/' {input.gp} | awk '/^set output/{{gsub(/.png/, ".large.png");}} /^set mouse/{{$0="# "$0;}} /^set style/{{$NF="0.5";}} 1' > {output.gp}
        cd $(dirname {output.gp})
        gnuplot $(basename {output.gp})
        ) &> {log}
        """
