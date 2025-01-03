rule helixer:
    input:
        "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa",
    output:
        "results/{asmname}/4.annotation/02.helixer/helixer.gff",
    log:
        "results/logs/4.annotation/helixer/{asmname}.log"
    benchmark:
        "results/benchmarks/4.annotation/helixer/{asmname}.txt"
    params:
        helixer_model = config["helixer_model"],
        subseqlen = config["helixer_max_gene_length"],
        species = lambda wildcards: get_species_name(wildcards),
    threads:
        lambda wildcards: len(config["reference_genomes"][get_reference_id(wildcards.asmname)]["chromosomes"]) + 1  #the number of chromosomes plus 1
    resources:
        helixer = 1
    container:
        "docker://gglyptodon/helixer-docker:helixer_v0.3.2_cuda_11.8.0-cudnn8"
    shell:
        "Helixer.py --fasta-path {input} --gff-output-path {output} --species {params.species} --subsequence-length {params.subseqlen} --model-filepath {params.helixer_model} &> {log}"
