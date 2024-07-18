include: "01.ntjoin.smk"
include: "02.renaming.smk"
include: "03.mummer.smk"

rule copy_assembly:
    input:
        "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa"
    output:
        protected("final_output/{asmname}.full.fa"),
    log:
        "results/logs/2.scaffolding/copy_assembly/{asmname}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/copy_assembly/{asmname}.txt"
    shell:
        "cp $(realpath {input}) {output} &> {log}"

def get_mummerplot_scaffolds(wildcards):
    filelist = []
    for asmname in get_all_accessions():
        reference = get_reference_id(asmname)
        filelist.append(f"results/{asmname}/2.scaffolding/03.mummer/{asmname}.vs.{reference}.plot.gp")
    return filelist

rule scaffold:
    input:
        expand("final_output/{asmname}.full.fa",
            asmname = get_all_accessions()),
        get_mummerplot_scaffolds,
    output:
        touch("results/scaffolding.done")
