rule hifiasm:
    input:
        hifi = lambda wildcards: config["reads"]["hifi"][wildcards.asmname],
    output:
        "results/{asmname}/1.assembly/01.hifiasm/{asmname}_hifi_only.p_ctg.gfa", #with `--primary -l0` option
    log:
        "results/logs/1.assembly/hifiasm/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/hifiasm/{asmname}.txt"
    threads:
        min(max(workflow.cores - 1, 1), 50)
    conda:
        "../../envs/hifiasm.yaml"
    shell:
        "hifiasm -t {threads} --primary -l0 -o $(echo {output} | rev | cut -d '.' -f 3- | rev) {input.hifi} 2> {log}"  # with `--primary -l0` option

rule hifiasm_with_ont:
    input:
        hifi = lambda wildcards: config["reads"]["hifi"][wildcards.asmname],
        ont = lambda wildcards: config["reads"]["ont"][wildcards.asmname],
    output:
        "results/{asmname}/1.assembly/01.hifiasm/{asmname}_hifi_and_ont.p_ctg.gfa", #with `--primary -l0` option
    log:
        "results/logs/1.assembly/hifiasm_with_ont/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/hifiasm_with_ont/{asmname}.txt"
    threads:
        min(max(workflow.cores - 1, 1), 50)
    conda:
        "../../envs/hifiasm.yaml"
    shell:
        "hifiasm -t {threads} --primary -l0 -o $(echo {output} | rev | cut -d '.' -f 3- | rev) --ul {input.ont} {input.hifi} 2> {log}"  # with `--primary -l0` option

rule gfa2fasta:
    input:
        lambda wildcards: branch("ont" in config["reads"] and wildcards.asmname in config["reads"]["ont"],  #with `--primary -l0` option
            then = "results/{asmname}/1.assembly/01.hifiasm/{asmname}_hifi_and_ont.p_ctg.gfa",
            otherwise = "results/{asmname}/1.assembly/01.hifiasm/{asmname}_hifi_only.p_ctg.gfa"),
    output:
        "results/{asmname}/1.assembly/01.hifiasm/{asmname}.fa"
    log:
        "results/logs/1.assembly/gfa2fasta/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/gfa2fasta/{asmname}.txt"
    shell:
        "awk '/^S/{{print \">\"$2; print $3;}}' {input} > {output} 2> {log}"
