rule sample_illumina:
    input:
        lambda wildcards: branch({wildcards.direction} == 1,
            then = get_illumina_1,
            otherwise = get_illumina_2),
    output:
        "results/{asmname}/5.quality_control/04.mapping/input/illumina_{direction}.fq.gz",
    log:
        "results/logs/5.quality_control/sample_illumina/{asmname}_{direction}.log"
    benchmark:
        "results/benchmarks/5.quality_control/sample_illumina/{asmname}_{direction}.txt"
    params:
        subset = 1000000
    conda:
        "../../envs/seqkit.yaml"
    shell:
        "seqkit head -n {params.subset} {input} -o {output} &> {log}"

rule bwa_index_genome:
    input:
        genome = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa",
    output:
        index1 = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa.0123",
        index2 = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa.amb",
        index3 = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa.ann",
        index4 = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa.bwt.2bit.64",
        index5 = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa.pac",
    log:
        "results/logs/5.quality_control/bwa_index_genome/{asmname}.log"
    benchmark:
        "results/benchmarks/5.quality_control/bwa_index_genome/{asmname}.txt"
    conda:
        "../../envs/mapping.yaml"
    shell:
        "bwa-mem2 index {input} &> {log}"

rule map_illumina:
    input:
        genome = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa",
        index1 = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa.0123",
        index2 = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa.amb",
        index3 = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa.ann",
        index4 = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa.bwt.2bit.64",
        index5 = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa.pac",
        forward = "results/{asmname}/5.quality_control/04.mapping/input/illumina_1.fq.gz",
        backward = "results/{asmname}/5.quality_control/04.mapping/input/illumina_2.fq.gz",
    output:
        "results/{asmname}/5.quality_control/04.mapping/output/illumina.sorted.bam",
    log:
        "results/logs/5.quality_control/map_illumina/{asmname}.log"
    benchmark:
        "results/benchmarks/5.quality_control/map_illumina/{asmname}.txt"
    params:
        readtag = lambda wildcards: f"'@RG\\tID:illumina\\tSM:{wildcards.asmname}'"
    threads:
        10
    conda:
        "../../envs/mapping.yaml"
    shell:
        "(bwa-mem2 mem -R {params.readtag} -t {threads} {input.genome} {input.forward} {input.backward} | samtools sort -@ $(({threads}-1)) -o {output} -) &> {log}"

rule index_bam:
    input:
        "results/{asmname}/5.quality_control/04.mapping/output/illumina.sorted.bam",
    output:
        "results/{asmname}/5.quality_control/04.mapping/output/illumina.sorted.bam.csi",
    log:
        "results/logs/5.quality_control/index_bam/{asmname}/illumina.log"
    benchmark:
        "results/benchmarks/5.quality_control/index_bam/{asmname}/illumina.txt"
    threads:
        10
    conda:
        "../../envs/mapping.yaml"
    shell:
        "samtools index -@ $(({threads}-1)) --csi {input} &> {log}"

rule flagstat_bam:
    input:
        bam = "results/{asmname}/5.quality_control/04.mapping/output/illumina.sorted.bam",
        index = "results/{asmname}/5.quality_control/04.mapping/output/illumina.sorted.bam.csi",
    output:
        "results/{asmname}/5.quality_control/04.mapping/output/illumina.sorted.flagstat.txt",
    log:
        "results/logs/5.quality_control/flagstat_bam/{asmname}/illumina.log"
    benchmark:
        "results/benchmarks/5.quality_control/flagstat_bam/{asmname}/illumina.txt"
    threads:
        10
    conda:
        "../../envs/mapping.yaml"
    shell:
        "samtools flagstat -@ $(({threads}-1)) {input.bam} > {output} 2> {log}"

def get_wgs_flagstat(wildcards):
    all_output = []
    if not SAMPLES[SAMPLES["accessionId"] == wildcards.asmname]["illumina_1"].isnull().values.item():
        all_output.append(f"results/{wildcards.asmname}/5.quality_control/04.mapping/output/illumina.sorted.flagstat.txt")
    return all_output

rule multiqc:
    input:
        get_wgs_flagstat
    output:
        report("results/{asmname}/5.quality_control/04.multiqc/multiqc_report.html",
            category="Mapping",
            caption="../../report/mapping.rst",
            labels={"type": "multiqc", "assembly": "{asmname}"}),
    log:
        "results/logs/5.quality_control/multiqc/{asmname}.log"
    benchmark:
        "results/benchmarks/5.quality_control/multiqc/{asmname}.txt"
    conda:
        "../../envs/multiqc.yaml"
    shell:
        "multiqc -f -o $(dirname {output}) $(dirname {input} | uniq) &> {log}"
