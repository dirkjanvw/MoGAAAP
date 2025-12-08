rule helixer_download_model:
    output:
        temporary(directory("helixer_models")),
    log:
        "results/logs/2.annotation/helixer_download_model.log"
    params:
        helixer_model = config["helixer_model"],
    container:
        "docker://gglyptodon/helixer-docker:helixer_v0.3.5_cuda_12.2.2-cudnn8"
    shell:
        """
        (
        mkdir {output}
        fetch_helixer_models.py --lineage {params.helixer_model} --custom-path {output}
        ) &> {log}
        """

rule helixer:
    input:
        genome = "final_output/{asmname}.full.fa",
        helixer_model = ancient(rules.helixer_download_model.output),
    output:
        "results/{asmname}/2.annotation/02.helixer/helixer.gff",
    log:
        "results/logs/2.annotation/helixer/{asmname}.log"
    benchmark:
        "results/benchmarks/2.annotation/helixer/{asmname}.txt"
    params:
        subseqlen = config["helixer_max_gene_length"],
        species = lambda wildcards: get_species_name(wildcards),
        helixer_model = config["helixer_model"],
    threads:
        lambda wildcards: len(config["reference_genomes"][get_reference_id(wildcards.asmname)]["chromosomes"]) + 1  #the number of chromosomes plus 1
    resources:
        helixer = 1
    container:
        "docker://gglyptodon/helixer-docker:helixer_v0.3.5_cuda_12.2.2-cudnn8"
    shell:
        "Helixer.py --fasta-path {input.genome} --gff-output-path {output} --species {params.species} --subsequence-length {params.subseqlen} --downloaded-model-path {input.helixer_model} --lineage {params.helixer_model} &> {log}"
