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
        filelist.append(f"results/{asmname}/2.scaffolding/03.mummer/{asmname}.vs.{reference}.plot.large.gp")
    return filelist

def get_hic_plots(wildcards):
    filelist = []
    for asmname in get_all_accessions():
        if has_hic(asmname):
            filelist.append(f"results/{asmname}/2.scaffolding/01.ntjoin/contact_map.pdf")
    return filelist

rule scaffold:
    input:
        expand("final_output/{asmname}.full.fa",
            asmname = get_all_accessions()),
        get_mummerplot_scaffolds,
        get_hic_plots,
    output:
        touch("results/scaffolding.done")
