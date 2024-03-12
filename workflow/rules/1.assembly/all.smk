include: "01.hifiasm.smk"
include: "02.contigs.smk"

rule assemble:
    input:
        expand("results/{asmname}/1.assembly/02.contigs/{asmname}.min{minlen}.sorted.renamed.fa",
                asmname=get_all_accessions(),
                minlen=config["min_contig_len"]
                ),
    output:
        touch("results/assembly.done")
