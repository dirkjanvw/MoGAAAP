rule renaming_sequences:
    input:
        reference = config["ref_genome"],
        scaffolds = expand("results/{{asmname}}/2.scaffolding/01.ntjoin/{{asmname}}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.assigned.scaffolds.fa",
                reference=config["ref_genome"],
                minlen=config["minlen"],
                k=config["k"],
                w=config["w"],
                ),
        delta = expand("results/{{asmname}}/2.scaffolding/02.mummer/{{asmname}}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.MUMmer.delta",
                reference=config["ref_genome"],
                minlen=config["minlen"],
                k=config["k"],
                w=config["w"],
                ),
    output:
        "results/{asmname}/2.scaffolding/03.renaming/{asmname}.fa"
    log:
        "results/logs/renaming_sequences/{asmname}.log"
    benchmark:
        "results/benchmarks/renaming_sequences/{asmname}.txt"
    shell:
        "??? &> {log}" #TODO: add script
