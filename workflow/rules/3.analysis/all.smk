include: "00.asm_bed.smk"
include: "01.blastdb.smk"
include: "02.blast_p.smk"
include: "03.blast_n.smk"
include: "04.blp2bed.smk"
include: "05.bln2bed.smk"
include: "06.bcovblp.smk"
include: "07.bcovbln.smk"
include: "08.circos.smk"
include: "09.separate_genome.smk"

rule link_separated_genome:
    input:
        nuclear = "results/{asmname}/3.analysis/09.separate_genome/{asmname}.nuclear.fa",
        organellar = expand("results/{{asmname}}/3.analysis/09.separate_genome/{{asmname}}.{organelle}.fa", organelle = config["organellar"]),
    output:
        nuclear = "final_output/{asmname}.nuclear.fa",
        organellar = expand("final_output/{{asmname}}.{organelle}.fa", organelle = config["organellar"]),
    log:
        "results/logs/3.analysis/link_separated_genome/{asmname}.log"
    benchmark:
        "results/benchmarks/3.analysis/link_separated_genome/{asmname}.txt"
    shell:
        """
        (
        ln -s $(realpath {input.nuclear}) {output.nuclear}
        ln -s $(realpath {input.organellar}) $(dirname {output.organellar} | uniq)
        ) &> {log}
        """

rule analyse:
    input:
        expand("results/{asmname}/3.analysis/08.circos/{asmname}.circos.html", asmname = get_all_accessions()), ### CIRCOS configuration ###
        expand("results/{asmname}/3.analysis/08.circos/{asmname}.circos.png", asmname = get_all_accessions()), ### CIRCOS PLOT ###
        expand("final_output/{asmname}.{section}.fa", section = [x for x in config["organellar"]] + ["nuclear"], asmname = get_all_accessions()), ### SEPARATED GENOMES ###
    output:
        touch("results/analysis.done")
