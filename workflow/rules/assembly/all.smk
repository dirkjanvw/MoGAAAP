include: "01.hifiasm.smk"

rule assemble:
    input:
        "results/{asmname}/assembly/01.hifiasm/{asmname}.min{minlen}.sorted.renamed.fa"
    output:
        touch("results/1.assembly/.done")
