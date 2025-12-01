include: "01.liftoff.smk"
include: "02.helixer.smk"
include: "03.combined.smk"
include: "04.asm_bed.smk"
include: "05.blastdb.smk"
include: "06.blast_p.smk"
include: "07.blast_n.smk"
include: "08.blp2bed.smk"
include: "09.bln2bed.smk"
include: "10.bcovblp.smk"
include: "11.bcovbln.smk"
include: "12.circos.smk"
include: "13.separate_genome.smk"
include: "14.telo.smk"

rule copy_annotation:
    input:
        full = lambda wildcards: get_annotation_location(wildcards.asmname),
        coding = "results/{asmname}/2.annotation/03.combined/{asmname}.coding.gff",
        clean = "results/{asmname}/2.annotation/03.combined/{asmname}.clean.gff",
    output:
        full = "final_output/{asmname}.full.gff",
        coding = "final_output/{asmname}.full.coding.gff",
        clean = "final_output/{asmname}.full.clean.gff",
    log:
        "results/logs/2.annotation/copy_annotation/{asmname}.log"
    benchmark:
        "results/benchmarks/2.annotation/copy_annotation/{asmname}.txt"
    shell:
        """
        (
        cp {input.full} {output.full}
        cp {input.coding} {output.coding}
        cp {input.clean} {output.clean}
        ) &> {log}
        """

rule copy_separated_genome:
    input:
        nuclear = "results/{asmname}/2.annotation/13.separate_genome/{asmname}.nuclear.fa",
        organellar = expand("results/{{asmname}}/2.annotation/13.separate_genome/{{asmname}}.{organelle}.fa", organelle = config["organellar"]),
    output:
        nuclear = "final_output/{asmname}.nuclear.fa",
        organellar = expand("final_output/{{asmname}}.{organelle}.fa", organelle = config["organellar"]),
    log:
        "results/logs/2.annotation/copy_separated_genome/{asmname}.log"
    benchmark:
        "results/benchmarks/2.annotation/copy_separated_genome/{asmname}.txt"
    shell:
        """
        (
        cp $(realpath {input.nuclear}) {output.nuclear}
        cp $(realpath {input.organellar}) $(dirname {output.organellar} | uniq)
        ) &> {log}
        """

def get_query_files(wildcards):
    query_files = []

    if "prot_queries" in config:
        for query_name in config["prot_queries"]:
            for asmname in get_all_accessions():
                query_files.append(f"results/{asmname}/2.annotation/06.blast_p/{query_name}.vs.{asmname}.html") ### BLASTP results ###

    if "nucl_queries" in config:
        for query_name in config["nucl_queries"]:
            for asmname in get_all_accessions():
                query_files.append(f"results/{asmname}/2.annotation/07.blast_n/{query_name}.vs.{asmname}.html") ### BLASTN results results for queries ###

    return query_files

rule annotate_genes:
    input:
        expand("final_output/{asmname}.full.gff", asmname=get_all_accessions()),
        expand("final_output/{asmname}.full.coding.gff", asmname=get_all_accessions()),
        expand("final_output/{asmname}.full.clean.gff", asmname=get_all_accessions()),
    output:
        touch("results/gene_annotation.done")

rule annotate_custom:
    input:
        get_query_files, ### optional BLAST results ###
        expand("results/{asmname}/2.annotation/07.blast_n/{query_name}.vs.{asmname}.html", query_name = config["organellar"], asmname = get_all_accessions()), ### BLASTN results for organellar ###
        expand("results/{asmname}/2.annotation/12.circos/{asmname}.circos.html", asmname = get_all_accessions()), ### CIRCOS configuration ###
        expand("results/{asmname}/2.annotation/12.circos/{asmname}.circos.png", asmname = get_all_accessions()), ### CIRCOS PLOT ###
        expand("results/{asmname}/2.annotation/14.telo/{asmname}.telo.html", asmname = get_all_accessions()), ### TELOMERE LOCATIONS ###
        expand("final_output/{asmname}.{section}.fa", section = [x for x in config["organellar"]] + ["nuclear"], asmname = get_all_accessions()), ### SEPARATED GENOMES ###
    output:
        touch("results/custom_annotation.done")

rule annotate:
    input:
        "results/gene_annotation.done",
        "results/custom_annotation.done",
    output:
        touch("results/annotation.done")
