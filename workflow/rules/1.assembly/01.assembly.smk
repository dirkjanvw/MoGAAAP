rule hifiasm:
    input:
        hifi = get_hifi,
    output:
        "results/{asmname}/1.assembly/01.hifiasm_hifi_only/{asmname}.bp.p_ctg.gfa",
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
        "results/{asmname}/1.assembly/01.hifiasm_hifi_and_ont/{asmname}.bp.p_ctg.gfa",
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
        "results/{asmname}/1.assembly/01.hifiasm_{ext}/{asmname}.bp.p_ctg.gfa",
    output:
        "results/{asmname}/1.assembly/01.hifiasm_{ext}/{asmname}.fa",
    log:
        "results/logs/1.assembly/hifiasm_{ext}_to_fasta/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/hifiasm_{ext}_to_fasta/{asmname}.txt"
    shell:
        "awk '/^S/{{print \">\"$2; print $3;}}' {input} > {output} 2> {log}"

rule flye:
    input:
        hifi = get_hifi,
    output:
        "results/{asmname}/1.assembly/01.flye_hifi_only/assembly.fasta",
    log:
        "results/logs/1.assembly/flye_hifi_only/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/flye_hifi_only/{asmname}.txt"
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
        "results/{asmname}/1.assembly/01.flye_hifi_and_ont/assembly.fasta",
    log:
        "results/logs/1.assembly/flye_hifi_and_ont/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/flye_hifi_and_ont/{asmname}.txt"
    threads:
        min(max(workflow.cores - 1, 1), 50)
    conda:
        "../../envs/flye.yaml"
    shell:
        "flye --pacbio-hifi {input.hifi} --nano-hq {input.ont} --out-dir $(dirname {output}) --threads {threads} &> {log}"

rule flye_rename_fasta:
    input:
        "results/{asmname}/1.assembly/01.flye_{ext}/assembly.fasta",
    output:
        "results/{asmname}/1.assembly/01.flye_{ext}/{asmname}.fa",
    log:
        "results/logs/1.assembly/flye_{ext}_rename_fasta/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/flye_{ext}_rename_fasta/{asmname}.txt"
    shell:
        "cp {input} {output} &> {log}"

rule verkko:
    input:
        hifi = get_hifi,
    output:
        "results/{asmname}/1.assembly/01.verkko_hifi_only/assembly.fasta",
    log:
        "results/logs/1.assembly/verkko_hifi_only/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/verkko_hifi_only/{asmname}.txt"
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
        "results/{asmname}/1.assembly/01.verkko_hifi_and_ont/assembly.fasta",
    log:
        "results/logs/1.assembly/verkko_hifi_and_ont/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/verkko_hifi_and_ont/{asmname}.txt"
    threads:
        min(max(workflow.cores - 1, 1), 50)
    conda:
        "../../envs/verkko.yaml"
    shell:
        "verkko --threads {threads} -d $(dirname {output}) --hifi {input.hifi} --nano {input.ont} &> {log}"

rule verkko_rename_fasta:
    input:
        "results/{asmname}/1.assembly/01.verkko_{ext}/assembly.fasta",
    output:
        "results/{asmname}/1.assembly/01.verkko_{ext}/{asmname}.fa",
    log:
        "results/logs/1.assembly/verkko_{ext}_rename_fasta/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/verkko_{ext}_rename_fasta/{asmname}.txt"
    shell:
        "cp {input} {output} &> {log}"
