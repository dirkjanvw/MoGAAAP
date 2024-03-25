rule sample_illumina:
    input:
        lambda wildcards: config["reads"]["illumina"][wildcards.asmname][wildcards.sample][wildcards.library][int(wildcards.direction)],
    output:
        "results/{asmname}/5.quality_control/04.mapping/input/illumina/{sample}/{library}_{direction}.fq.gz",
    log:
        "results/logs/5.quality_control/sample_illumina/{asmname}/{sample}/{library}_{direction}.log"
    benchmark:
        "results/benchmarks/5.quality_control/sample_illumina/{asmname}/{sample}/{library}_{direction}.txt"
    params:
        subset = 1000000
    conda:
        "../../envs/seqkit.yaml"
    shell:
        "seqkit head -n {params.subset} {input} -o {output} &> {log}"

rule map_illumina:
    input:
        genome = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa",
        index1 = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa.0123",
        index2 = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa.amb",
        index3 = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa.ann",
        index4 = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa.bwt.2bit.64",
        index5 = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa.pac",
        forward = "results/{asmname}/5.quality_control/04.mapping/input/illumina/{sample}/{library}_1.fq.gz",
        backward = "results/{asmname}/5.quality_control/04.mapping/input/illumina/{sample}/{library}_2.fq.gz",
    output:
        "results/{asmname}/5.quality_control/04.mapping/output/illumina/{sample}/{library}.sorted.bam",
    log:
        "results/logs/5.quality_control/map_illumina/{asmname}/{sample}/{library}.log"
    benchmark:
        "results/benchmarks/5.quality_control/map_illumina/{asmname}/{sample}/{library}.txt"
    params:
        readtag = lambda wildcards: f"'@RG\\tID:{wildcards.library}\\tSM:{wildcards.library}'"
    threads:
        10
    conda:
        "../../envs/mapping.yaml"
    shell:
        "(bwa-mem2 mem -R {params.readtag} -t {threads} {input.genome} {input.forward} {input.backward} | samtools sort -@ $(({threads}-1)) -o {output} -) &> {log}"

rule sample_tgs:
    input:
        lambda wildcards: config["reads"][wildcards.tgstype][wildcards.asmname][wildcards.sample]
    output:
        "results/{asmname}/5.quality_control/04.mapping/input/{tgstype}/{sample}/{sample}.fq.gz"
    log:
        "results/logs/5.quality_control/sample_tgs/{asmname}/{tgstype}/{sample}.log"
    benchmark:
        "results/benchmarks/5.quality_control/sample_tgs/{asmname}/{tgstype}/{sample}.txt"
    params:
        subset = 10000
    conda:
        "../../envs/seqkit.yaml"
    shell:
        "seqkit head -n {params.subset} {input} -o {output} &> {log}"

rule map_tgs:
    input:
        genome = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa",
        reads = "results/{asmname}/5.quality_control/04.mapping/input/{tgstype}/{sample}/{sample}.fq.gz",
    output:
        "results/{asmname}/5.quality_control/04.mapping/output/{tgstype}/{sample}/{sample}.sorted.bam",
    log:
        "results/logs/5.quality_control/map_tgs/{asmname}/{tgstype}/{sample}.log"
    benchmark:
        "results/benchmarks/5.quality_control/map_tgs/{asmname}/{tgstype}/{sample}.txt"
    params:
        readtag = lambda wildcards: f"'@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}'"
    threads:
        10
    conda:
        "../../envs/mapping.yaml"
    shell:
        "(minimap2 -ax map-{wildcards.tgstype} -t {threads} -R {params.readtag} {input.genome} {input.reads} | samtools sort -@ $(({threads}-1)) -o {output} -) &> {log}"

rule index_bam:
    input:
        "results/{asmname}/5.quality_control/04.mapping/output/{wgstype}/{sample}/{library}.sorted.bam",
    output:
        "results/{asmname}/5.quality_control/04.mapping/output/{wgstype}/{sample}/{library}.sorted.bam.csi",
    log:
        "results/logs/5.quality_control/index_bam/{asmname}/{wgstype}/{sample}/{library}.log"
    benchmark:
        "results/benchmarks/5.quality_control/index_bam/{asmname}/{wgstype}/{sample}/{library}.txt"
    threads:
        10
    conda:
        "../../envs/mapping.yaml"
    shell:
        "samtools index -@ $(({threads}-1)) --csi {input} &> {log}"

rule flagstat_bam:
    input:
        bam = "results/{asmname}/5.quality_control/04.mapping/output/{wgstype}/{sample}/{library}.sorted.bam",
        index = "results/{asmname}/5.quality_control/04.mapping/output/{wgstype}/{sample}/{library}.sorted.bam.csi",
    output:
        "results/{asmname}/5.quality_control/04.mapping/output/{wgstype}/{sample}/{library}.sorted.flagstat.txt",
    log:
        "results/logs/5.quality_control/flagstat_bam/{asmname}/{wgstype}/{sample}/{library}.log"
    benchmark:
        "results/benchmarks/5.quality_control/flagstat_bam/{asmname}/{wgstype}/{sample}/{library}.txt"
    threads:
        10
    conda:
        "../../envs/mapping.yaml"
    shell:
        "samtools flagstat -@ $(({threads}-1)) {input.bam} > {output} 2> {log}"

def get_wgs_flagstat(wildcards):
    all_output = []
    for sample in config["reads"]["hifi"][wildcards.asmname]:
        all_output.append(f"results/{wildcards.asmname}/5.quality_control/04.mapping/output/hifi/{sample}/{sample}.sorted.flagstat.txt")
    if "ont" in config["reads"]:
        if wildcards.asmname in config["reads"]["ont"]:
            for sample in config["reads"]["ont"][wildcards.asmname]:
                all_output.append(f"results/{wildcards.asmname}/5.quality_control/04.mapping/output/ont/{sample}/{sample}.sorted.flagstat.txt")
    if "illumina" in config["reads"]:
        if wildcards.asmname in config["reads"]["illumina"]:
            for sample in config["reads"]["illumina"][wildcards.asmname]:
                for library in config["reads"]["illumina"][wildcards.asmname][sample]:
                    all_output.append(f"results/{wildcards.asmname}/5.quality_control/04.mapping/output/illumina/{sample}/{library}.sorted.flagstat.txt")
    return all_output

rule multiqc:
    input:
        get_wgs_flagstat
    output:
        report("results/{asmname}/5.quality_control/04.multiqc/multiqc_report.html", category="Mapping", labels={"type": "multiqc", "assembly": "{asmname}"}),
    log:
        "results/logs/5.quality_control/multiqc/{asmname}.log"
    benchmark:
        "results/benchmarks/5.quality_control/multiqc/{asmname}.txt"
    conda:
        "../../envs/multiqc.yaml"
    shell:
        "multiqc -f -o $(dirname {output}) $(dirname {input} | uniq) &> {log}"
