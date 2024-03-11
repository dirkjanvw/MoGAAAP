include: "01.ntjoin.smk"
include: "02.renaming.smk"
include: "03.mummer.smk"

rule scaffold:
    input:
        expand("results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa",
            asmname = get_all_accessions()),
        expand("results/{asmname}/2.scaffolding/03.mummer/{asmname}.vs.{reference}.plot.gp",
            asmname = get_all_accessions(),
            reference = config["ref_genome"]),
    output:
        touch("results/scaffolding.done")
