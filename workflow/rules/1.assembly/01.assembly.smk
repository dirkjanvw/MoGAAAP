rule hifiasm:
    input:
        hifi = get_hifi,
    output:
        hap1 = "results/{asmname}/1.assembly/01.hifiasm_hifi_only/{asmname}.bp.hap1.p_ctg.gfa",
        hap2 = "results/{asmname}/1.assembly/01.hifiasm_hifi_only/{asmname}.bp.hap2.p_ctg.gfa",
        consensus = "results/{asmname}/1.assembly/01.hifiasm_hifi_only/{asmname}.bp.p_ctg.gfa",
    log:
        "results/logs/1.assembly/hifiasm/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/hifiasm/{asmname}.txt"
    threads:
        min(max(workflow.cores - 1, 1), 50)
    conda:
        "../../envs/hifiasm.yaml"
    shell:
        "hifiasm -t {threads} -o $(echo {output.consensus} | rev | cut -d '.' -f 4- | rev) {input.hifi} 2> {log}"

rule hifiasm_with_ont:
    input:
        hifi = get_hifi,
        ont = get_ont,
    output:
        hap1 = "results/{asmname}/1.assembly/01.hifiasm_hifi_and_ont/{asmname}.bp.hap1.p_ctg.gfa",
        hap2 = "results/{asmname}/1.assembly/01.hifiasm_hifi_and_ont/{asmname}.bp.hap2.p_ctg.gfa",
        consensus = "results/{asmname}/1.assembly/01.hifiasm_hifi_and_ont/{asmname}.bp.p_ctg.gfa",
    log:
        "results/logs/1.assembly/hifiasm_with_ont/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/hifiasm_with_ont/{asmname}.txt"
    threads:
        min(max(workflow.cores - 1, 1), 50)
    conda:
        "../../envs/hifiasm.yaml"
    shell:
        "hifiasm -t {threads} -o $(echo {output.consensus} | rev | cut -d '.' -f 4- | rev) --ul $(echo {input.ont} | sed 's/ /,/g') {input.hifi} 2> {log}"

def get_hifiasm_output(wildcards):
    if get_haplotypes(wildcards) == 1:
        return "results/{asmname}/1.assembly/01.hifiasm_{ext}/{asmname}.bp.p_ctg.gfa"
    else:
        haplotype = int(wildcards.asmname[-1])
        asmname = get_clean_accession_id(wildcards.asmname)
        return f"results/{asmname}/1.assembly/01.hifiasm_{wildcards.ext}/{asmname}.bp.hap{haplotype}.p_ctg.gfa"

rule hifiasm_to_fasta:
    input:
        get_hifiasm_output,
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

rule other_assembler_assembly_minimap2:
    input:
        hifi=get_hifi,
        assembly="results/{asmname}/1.assembly/01.{assembler}_{ext}/assembly.fasta",
    output:
        "results/{asmname}/1.assembly/01.{assembler}_{ext}/assembly_hapdup/hifi_vs_assembly.bam",
    log:
        "results/logs/1.assembly/{assembler}_{ext}_assembly_minimap2/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/{assembler}_{ext}_assembly_minimap2/{asmname}.txt"
    threads:
        min(max(workflow.cores - 1,1),50)
    conda:
        "../../envs/minimap2.yaml"
    shell:
        "(minimap2 -ax map-hifi -t {threads} {input.assembly} {input.hifi} | samtools sort -@ $(({threads}-1)) > {output}) 2> {log}"

rule other_assembler_assembly_minimap2_index:
    input:
        bam="results/{asmname}/1.assembly/01.{assembler}_{ext}/assembly_hapdup/hifi_vs_assembly.bam",
    output:
        "results/{asmname}/1.assembly/01.{assembler}_{ext}/assembly_hapdup/hifi_vs_assembly.bam.bai",
    log:
        "results/logs/1.assembly/{assembler}_{ext}_assembly_minimap2_index/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/{assembler}_{ext}_assembly_minimap2_index/{asmname}.txt"
    threads:
        min(max(workflow.cores - 1,1),50)
    conda:
        "../../envs/samtools.yaml"
    shell:
        "samtools index -@ $(({threads}-1)) {input} &> {log}"

rule other_assembler_hapdup:
    input:
        bam="results/{asmname}/1.assembly/01.{assembler}_{ext}/assembly_hapdup/hifi_vs_assembly.bam",
        bai="results/{asmname}/1.assembly/01.{assembler}_{ext}/assembly_hapdup/hifi_vs_assembly.bam.bai",
        assembly="results/{asmname}/1.assembly/01.{assembler}_{ext}/assembly.fasta",
    output:
        "results/{asmname}/1.assembly/01.{assembler}_{ext}/assembly_hapdup/hapdup_dual_1.fasta",
        "results/{asmname}/1.assembly/01.{assembler}_{ext}/assembly_hapdup/hapdup_dual_2.fasta",
        "results/{asmname}/1.assembly/01.{assembler}_{ext}/assembly_hapdup/phased_blocks_hp1.bed",
        "results/{asmname}/1.assembly/01.{assembler}_{ext}/assembly_hapdup/phased_blocks_hp2.bed",
        "results/{asmname}/1.assembly/01.{assembler}_{ext}/assembly_hapdup/hapdup_phased_1.fasta",
        "results/{asmname}/1.assembly/01.{assembler}_{ext}/assembly_hapdup/hapdup_phased_2.fasta",
    log:
        "results/logs/1.assembly/{assembler}_{ext}_hapdup/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/{assembler}_{ext}_hapdup/{asmname}.txt"
    threads:
        min(max(workflow.cores - 1,1),50)
    container:
        "docker://mkolmogo/hapdup:0.12"
    shell:
        "hapdup --threads {threads} --bam {input.bam} --rtype hifi --assembly {input.assembly} --out-dir $(dirname {output[0]}) &> {log}"

def get_other_assembler_output(wildcards):
    if get_haplotypes(wildcards) == 1:
        return "results/{asmname}/1.assembly/01.{assembler}_{ext}/assembly.fasta"
    else:
        haplotype = int(wildcards.asmname[-1])
        asmname = get_clean_accession_id(wildcards.asmname)
        return f"results/{asmname}/1.assembly/01.{{assembler}}_{{ext}}/assembly_hapdup/hapdup_dual_{haplotype}.fasta"  #using dual output here instead of phased, since we would like to keep contiguity (so phase switching may occur)

rule other_assembler_rename_fasta:
    input:
        get_other_assembler_output,
    output:
        "results/{asmname}/1.assembly/01.{assembler}_{ext}/{asmname}.fa",
    log:
        "results/logs/1.assembly/{assembler}_{ext}_rename_fasta/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/{assembler}_{ext}_rename_fasta/{asmname}.txt"
    shell:
        "cp {input} {output} &> {log}"
