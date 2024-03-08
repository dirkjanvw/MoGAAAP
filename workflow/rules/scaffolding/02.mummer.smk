rule mummer:
    input:
        reference = lambda wildcards: config["ref_genome"][wildcards.reference],
        scaffolds = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.assigned.scaffolds.fa",
    output:
        "results/{asmname}/2.scaffolding/02.mummer/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.MUMmer.delta",
    log:
        "results/logs/mummer/{reference}/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.log"
    benchmark:
        "results/benchmarks/mummer/{reference}/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.txt"
    threads:
        24
    conda:
        "../../envs/mummer.yaml"
    shell:
        "nucmer -t {threads} -l 1000 -g 1000 --prefix=$(echo {output} | rev | cut -d '.' -f 2- | rev) {input.reference} {input.scaffolds} &> {log}"

rule dotplot:
    input:
        "results/{asmname}/2.scaffolding/02.mummer/{alignment}.MUMmer.delta",
    output:
        "results/{asmname}/2.scaffolding/02.mummer/{alignment}.MUMmer.plot.gp",
    log:
        "results/logs/dotplot/{asmname}/{alignment}.log"
    benchmark:
        "results/benchmarks/dotplot/{asmname}/{alignment}.txt"
    conda:
        "../../envs/mummer.yaml"
    shell:
        """
        (
        ln -s $(realpath {input}) $(dirname {output})/
        cd $(dirname {output})
        mummerplot --filter --png --large --prefix={wildcards.alignment}.MUMmer.plot --title {wildcards.asmname} $(basename {input})
        ) &> {log}
        """