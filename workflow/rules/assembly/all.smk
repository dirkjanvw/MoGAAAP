include: "01.hifiasm.smk"
include: "02.contigs.smk"

rule assemble:
    input:
        "results/{asmname}/1.assembly/01.hifiasm/{asmname}.min{minlen}.sorted.renamed.fa"
    output:
        touch("results/1.assembly/.done")
