include: "00.asm_bed.smk"
include: "01.blastdb.smk"
include: "02.blast_p.smk"
include: "03.blast_n.smk"
include: "04.blp2bed.smk"
include: "05.bln2bed.smk"
include: "06.bcovblp.smk"
include: "07.bcovbln.smk"

rule analyse:
    input:
        expand("results/{asmname}/3.analysis/08.circos/{asmname}.circos.png", asmname = get_all_accessions()), ### CIRCOS PLOT ###
    output:
        touch("results/3.analysis/.done")
