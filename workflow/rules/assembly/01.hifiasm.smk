rule hifiasm:
    input:
        hifi = lambda wildcards: config["reads"]["hifi"][wildcards.asmname],
    output:
        "results/assembly/{asmname}/01.hifiasm/{asmname}_hifi_only.bp.p_ctg.gfa", #without `--primary -l0` option
        "results/assembly/{asmname}/01.hifiasm/{asmname}_hifi_only.p_ctg.gfa", #with `--primary -l0` option
    log:
        "results/logs/01.hifiasm/hifiasm/{asmname}.log"
    benchmark:
        "results/benchmarks/01.hifiasm/hifiasm/{asmname}.txt"
    threads:
        min(workflow.cores - 1, 50)
    conda:
        "../../envs/hifiasm.yaml"
    shell:
        "hifiasm -t {threads} -o $(echo {output} | rev | cut -d '.' -f 4- | rev) {input.hifi} 2> {log}"  # without `--primary -l0` option
        "hifiasm -t {threads} --primary -l0 -o $(echo {output} | rev | cut -d '.' -f 4- | rev) {input.hifi} 2> {log}"  # with `--primary -l0` option

rule hifiasm_with_ont:
    input:
        hifi = lambda wildcards: config["reads"]["hifi"][wildcards.asmname],
        ont = lambda wildcards: config["reads"]["ont"][wildcards.asmname],
    output:
        "results/assembly/{asmname}/01.hifiasm/{asmname}_hifi_and_ont.bp.p_ctg.gfa", #without `--primary -l0` option
        "results/assembly/{asmname}/01.hifiasm/{asmname}_hifi_and_ont.p_ctg.gfa", #with `--primary -l0` option
    log:
        "results/logs/01.hifiasm/hifiasm_with_ont/{asmname}.log"
    benchmark:
        "results/benchmarks/01.hifiasm/hifiasm_with_ont/{asmname}.txt"
    threads:
        min(workflow.cores - 1, 50)
    conda:
        "../../envs/hifiasm.yaml"
    shell:
        "hifiasm -t {threads} -o $(echo {output} | rev | cut -d '.' -f 4- | rev) --ul {input.ont} {input.hifi} 2> {log}"  # without `--primary -l0` option
        "hifiasm -t {threads} --primary -l0 -o $(echo {output} | rev | cut -d '.' -f 4- | rev) --ul {input.ont} {input.hifi}  2> {log}"  # with `--primary -l0` option

rule gfa2fasta:
    input:
        lambda wildcards: branch(wildcards.asmname in config["reads"]["ont"],  #without `--primary -l0` option
            then = "results/assembly/{asmname}/01.hifiasm/{asmname}_hifi_and_ont.bp.p_ctg.gfa",
            otherwise = "results/assembly/{asmname}/01.hifiasm/{asmname}_hifi_only.bp.p_ctg.gfa"),
        lambda wildcards: branch(wildcards.asmname in config["reads"]["ont"],  #with `--primary -l0` option
            then = "results/assembly/{asmname}/01.hifiasm/{asmname}_hifi_and_ont.p_ctg.gfa",
            otherwise = "results/assembly/{asmname}/01.hifiasm/{asmname}_hifi_only.p_ctg.gfa"),
    output:
        "results/{asmname}/assembly/01.hifiasm/{asmname}.fa"
    log:
        "results/logs/03.contigs_all/gfa2fasta/{asmname}.log"
    benchmark:
        "results/benchmarks/03.contigs_all/gfa2fasta/{asmname}.txt"
    shell:
        "awk '/^S/{{print \">\"$2; print $3;}}' {input} > {output} 2> {log}"
