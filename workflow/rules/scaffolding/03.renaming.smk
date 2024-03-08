rule renaming_sequences:
    input:
        reference = lambda wildcards: config["ref_genome"][wildcards.reference],
        scaffolds = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.assigned.scaffolds.fa",
        delta = "results/{asmname}/2.scaffolding/02.mummer/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.MUMmer.delta",
    output:
        "results/{asmname}/2.scaffolding/03.renaming/{asmname}.fa"
    log:
        "results/logs/renaming_sequences/{asmname}.log"
    benchmark:
        "results/benchmarks/renaming_sequences/{asmname}.txt"
    shell:
        "??? &> {log}" #TODO: add script
