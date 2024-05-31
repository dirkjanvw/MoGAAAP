rule hifiasm:
    input:
        hifi = get_hifi,
    output:
        "results/{asmname}/1.assembly/01.hifiasm/{asmname}_hifi_only.bp.p_ctg.gfa",
    log:
        "results/logs/1.assembly/hifiasm/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/hifiasm/{asmname}.txt"
    threads:
        min(max(workflow.cores - 1, 1), 50)
    conda:
        "../../envs/hifiasm.yaml"
    shell:
        "hifiasm -t {threads} -o $(echo {output} | rev | cut -d '.' -f 4- | rev) {input.hifi} 2> {log}"

rule hifiasm_with_ont:
    input:
        hifi = get_hifi,
        ont = get_ont,
    output:
        "results/{asmname}/1.assembly/01.hifiasm/{asmname}_hifi_and_ont.bp.p_ctg.gfa",
    log:
        "results/logs/1.assembly/hifiasm_with_ont/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/hifiasm_with_ont/{asmname}.txt"
    threads:
        min(max(workflow.cores - 1, 1), 50)
    conda:
        "../../envs/hifiasm.yaml"
    shell:
        "hifiasm -t {threads} -o $(echo {output} | rev | cut -d '.' -f 4- | rev) --ul {input.ont} {input.hifi} 2> {log}"

rule gfa2fasta:
    input:
        lambda wildcards: branch(SAMPLES[SAMPLES["accessionId"] == wildcards.asmname]["ont"].isnull().values.item(),  #check if ont is null
            then = "results/{asmname}/1.assembly/01.hifiasm/{asmname}_hifi_only.bp.p_ctg.gfa",
            otherwise = "results/{asmname}/1.assembly/01.hifiasm/{asmname}_hifi_and_ont.bp.p_ctg.gfa")
    output:
        "results/{asmname}/1.assembly/01.hifiasm/{asmname}.fa"
    log:
        "results/logs/1.assembly/gfa2fasta/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/gfa2fasta/{asmname}.txt"
    shell:
        "awk '/^S/{{print \">\"$2; print $3;}}' {input} > {output} 2> {log}"
