include: "01.ntjoin.smk"
include: "02.mummer.smk"
include: "03.renaming.smk"

rule scaffold:
    input:
        expand("results/{asmname}/2.scaffolding/03.renaming/{asmname}.fa", asmname = get_all_accessions()),
    output:
        touch("results/2.scaffolding/.done")
