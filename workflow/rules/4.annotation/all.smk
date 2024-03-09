include: "01.liftoff.smk"
include: "02.helixer.smk"
include: "03.combined.smk"

rule annotate:
    input:
        expand("results/{asmname}/4.annotation/03.combined/{asmname}.gff", asmname=get_all_accessions()),
    output:
        touch("results/4.annotation/.done")
