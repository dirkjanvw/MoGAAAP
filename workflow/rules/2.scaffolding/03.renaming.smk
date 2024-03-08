rule renaming_sequences:
    input:
        reference = lambda wildcards: config["ref_genome"][wildcards.reference],
        scaffolds = lambda wildcards: expand("results/{{asmname}}/2.scaffolding/01.ntjoin/{{asmname}}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.assigned.scaffolds.fa",
                reference=config["ref_genome"][wildcards.reference],
                minlen=config["min_contig_len"],
                k=config["ntjoin_k"],
                w=config["ntjoin_w"],
                ),
        delta = lambda wildcards: expand("results/{{asmname}}/2.scaffolding/02.mummer/{{asmname}}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.MUMmer.delta",
                reference=config["ref_genome"][wildcards.reference],
                minlen=config["min_contig_len"],
                k=config["ntjoin_k"],
                w=config["ntjoin_w"],
                ),
    output:
        "results/{asmname}/2.scaffolding/03.renaming/{asmname}.fa"
    log:
        "results/logs/2.scaffolding/renaming_sequences/{asmname}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/renaming_sequences/{asmname}.txt"
    shell:
        "??? &> {log}" #TODO: add script

rule index_sequences:
    input:
        "results/{asmname}/2.scaffolding/03.renaming/{asmname}.fa"
    output:
        "results/{asmname}/2.scaffolding/03.renaming/{asmname}.fa.fai"
    log:
        "results/logs/2.scaffolding/index_sequences/{asmname}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/index_sequences/{asmname}.txt"
    conda:
        "../../envs/samtools.yaml"
    shell:
        "samtools faidx {input} &> {log}"
