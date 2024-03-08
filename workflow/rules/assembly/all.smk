include: "01.hifiasm.smk"

rule assemble:
    input:
        ""
    output:
        touch("results/1.assembly/.done")
