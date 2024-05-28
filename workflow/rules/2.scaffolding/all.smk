include: "01.ntjoin.smk"
include: "02.renaming.smk"
include: "03.mummer.smk"

def get_mummerplot_scaffolds(wildcards):
    filelist = []
    for asmname in get_all_accessions():
        reference = get_reference_id(asmname)
        filelist.append(f"results/{asmname}/2.scaffolding/03.mummer/{asmname}.vs.{reference}.plot.gp")
    return filelist

rule scaffold:
    input:
        expand("results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa",
            asmname = get_all_accessions()),
        get_mummerplot_scaffolds,
    output:
        touch("results/scaffolding.done")
