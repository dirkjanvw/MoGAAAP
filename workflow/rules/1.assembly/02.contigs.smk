rule filter_contigs:
    input:
        lambda wildcards: branch(SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(wildcards.asmname)]["ont"].isnull().values.item(), #check if ont is null
                then=f"results/{{asmname}}/1.assembly/01.{config['assembler']}_hifi_only/{{asmname}}.fa",
                otherwise=f"results/{{asmname}}/1.assembly/01.{config['assembler']}_hifi_and_ont/{{asmname}}.fa"),
    output:
        "results/{asmname}/1.assembly/02.contigs/{asmname}.min{minlen}.fa"
    log:
        "results/logs/1.assembly/filter_contigs/{asmname}.{minlen}.log"
    benchmark:
        "results/benchmarks/1.assembly/filter_contigs/{asmname}.{minlen}.txt"
    conda:
        "../../envs/seqkit.yaml"
    shell:
        "seqkit seq -m {wildcards.minlen} {input} > {output} 2> {log}"

rule sort_contigs:
    input:
        "results/{asmname}/1.assembly/02.contigs/{asmname}.min{minlen}.fa"
    output:
        "results/{asmname}/1.assembly/02.contigs/{asmname}.min{minlen}.sorted.fa"
    log:
        "results/logs/1.assembly/sort_contigs/{asmname}.min{minlen}.log"
    benchmark:
        "results/benchmarks/1.assembly/sort_contigs/{asmname}.min{minlen}.txt"
    conda:
        "../../envs/seqkit.yaml"
    shell:
        "seqkit sort -rl {input} > {output} 2> {log}"

rule add_prefix:
    input:
        "results/{asmname}/1.assembly/02.contigs/{asmname}.min{minlen}.sorted.fa"
    output:
        "results/{asmname}/1.assembly/02.contigs/{asmname}.min{minlen}.sorted.renamed.fa"
    log:
        "results/logs/1.assembly/add_prefix/{asmname}.min{minlen}.log"
    benchmark:
        "results/benchmarks/1.assembly/add_prefix/{asmname}.min{minlen}.txt"
    params:
        prefix=lambda wildcards: f"{wildcards.asmname}_C"
    conda:
        "../../envs/bioawk.yaml"
    shell:
        "bioawk -c fastx '{{ print \">{params.prefix}\" ++i; print $seq }}' {input} > {output} 2> {log}"

rule index_contigs:
    input:
        "results/{asmname}/1.assembly/02.contigs/{asmname}.min{minlen}.sorted.renamed.fa"
    output:
        "results/{asmname}/1.assembly/02.contigs/{asmname}.min{minlen}.sorted.renamed.fa.fai"
    log:
        "results/logs/1.assembly/index_contigs/{asmname}.min{minlen}.log"
    benchmark:
        "results/benchmarks/1.assembly/index_contigs/{asmname}.min{minlen}.txt"
    conda:
        "../../envs/samtools.yaml"
    shell:
        "samtools faidx {input} &> {log}"