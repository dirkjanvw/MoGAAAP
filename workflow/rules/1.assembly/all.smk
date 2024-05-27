include: "01.hifiasm.smk"
include: "02.contigs.smk"
include: "03.mummer.smk"

rule assemble:
    input:
        expand("results/{asmname}/1.assembly/02.contigs/{asmname}.min{minlen}.sorted.renamed.fa",
            asmname=get_all_accessions(),
            minlen=config["min_contig_len"]
        ),
        expand("results/{asmname}/1.assembly/03.mummer/{asmname}.min{minlen}.vs.reference.plot.gp",
            asmname = get_all_accessions(),
            minlen=config["min_contig_len"]
        ),
    output:
        touch("results/assembly.done")
