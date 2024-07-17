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
        "hifiasm -t {threads} -o $(echo {output} | rev | cut -d '.' -f 4- | rev) --ul $(echo {input.ont} | sed 's/ /,/g') {input.hifi} 2> {log}"

rule hifiasm_to_fasta:
    input:
        lambda wildcards: branch(SAMPLES[SAMPLES["accessionId"] == wildcards.asmname]["ont"].isnull().values.item(),  #check if ont is null
            then = "results/{asmname}/1.assembly/01.hifiasm/{asmname}_hifi_only.bp.p_ctg.gfa",
            otherwise = "results/{asmname}/1.assembly/01.hifiasm/{asmname}_hifi_and_ont.bp.p_ctg.gfa")
    output:
        "results/{asmname}/1.assembly/01.hifiasm/{asmname}.fa"
    log:
        "results/logs/1.assembly/hifiasm_to_fasta/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/hifiasm_to_fasta/{asmname}.txt"
    shell:
        "awk '/^S/{{print \">\"$2; print $3;}}' {input} > {output} 2> {log}"

rule flye:
    input:
        hifi = get_hifi,
    output:
        "results/{asmname}/1.assembly/01.flye/{asmname}_flye.fasta",
    log:
        "results/logs/1.assembly/flye/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/flye/{asmname}.txt"
    threads:
        min(max(workflow.cores - 1, 1), 50)
    conda:
        "../../envs/flye.yaml"
    shell:
        "flye --pacbio-hifi {input.hifi} --out-dir $(dirname {output}) --threads {threads} &> {log}"

rule flye_with_ont:
    input:
        hifi = get_hifi,
        ont = get_ont,
    output:
        "results/{asmname}/1.assembly/01.flye/{asmname}_flye_with_ont.fasta",
    log:
        "results/logs/1.assembly/flye_with_ont/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/flye_with_ont/{asmname}.txt"
    threads:
        min(max(workflow.cores - 1, 1), 50)
    conda:
        "../../envs/flye.yaml"
    shell:
        "flye --pacbio-hifi {input.hifi} --nano-hq {input.ont} --out-dir $(dirname {output}) --threads {threads} &> {log}"

rule flye_to_fasta:
    input:
        lambda wildcards: branch(SAMPLES[SAMPLES["accessionId"] == wildcards.asmname]["ont"].isnull().values.item(),  #check if ont is null
            then = "results/{asmname}/1.assembly/01.flye/{asmname}_flye.fasta",
            otherwise = "results/{asmname}/1.assembly/01.flye/{asmname}_flye_with_ont.fasta")
    output:
        "results/{asmname}/1.assembly/01.flye/{asmname}.fa"
    log:
        "results/logs/1.assembly/flye_to_fasta/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/flye_to_fasta/{asmname}.txt"
    shell:
        "awk '/^>/{print;next}{print toupper($0)}' {input} > {output} 2> {log}"

rule verkko:
    input:
        hifi = get_hifi,
    output:
        "results/{asmname}/1.assembly/01.verkko/{asmname}_verkko.fasta",
    log:
        "results/logs/1.assembly/verkko/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/verkko/{asmname}.txt"
    threads:
        min(max(workflow.cores - 1, 1), 50)
    conda:
        "../../envs/verkko.yaml"
    shell:
        "verkko --threads {threads} -d $(dirname {output}) --hifi {input.hifi} &> {log}"

rule verkko_with_ont:
    input:
        hifi = get_hifi,
        ont = get_ont,
    output:
        "results/{asmname}/1.assembly/01.verkko/{asmname}_verkko_with_ont.fasta",
    log:
        "results/logs/1.assembly/verkko_with_ont/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/verkko_with_ont/{asmname}.txt"
    threads:
        min(max(workflow.cores - 1, 1), 50)
    conda:
        "../../envs/verkko.yaml"
    shell:
        "verkko --threads {threads} -d $(dirname {output}) --hifi {input.hifi} --nano {input.ont} &> {log}"

rule verkko_to_fasta:
    input:
        lambda wildcards: branch(SAMPLES[SAMPLES["accessionId"] == wildcards.asmname]["ont"].isnull().values.item(),  #check if ont is null
            then = "results/{asmname}/1.assembly/01.verkko/{asmname}_verkko.fasta",
            otherwise = "results/{asmname}/1.assembly/01.verkko/{asmname}_verkko_with_ont.fasta")
    output:
        "results/{asmname}/1.assembly/01.verkko/{asmname}.fa"
    log:
        "results/logs/1.assembly/verkko_to_fasta/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/verkko_to_fasta/{asmname}.txt"
    shell:
        "awk '/^>/{print;next}{print toupper($0)}' {input} > {output} 2> {log}"
