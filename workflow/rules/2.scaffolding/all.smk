include: "01.scaffolding.smk"
include: "02.renaming.smk"
include: "03.mummer.smk"

rule copy_assembly:
    input:
        lambda wildcards:
            "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa"
            if not has_assembly_location(wildcards.asmname)
            else get_assembly_location(wildcards.asmname),
    output:
        "final_output/{asmname}.full.fa",
    log:
        "results/logs/2.scaffolding/copy_assembly/{asmname}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/copy_assembly/{asmname}.txt"
    shell:
        "cp $(realpath {input}) {output} &> {log}"

rule index_scaffolds:
    input:
        "final_output/{asmname}.full.fa"
    output:
        "final_output/{asmname}.full.fa.fai"
    log:
        "results/logs/2.scaffolding/index_scaffolds/{asmname}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/index_scaffolds/{asmname}.txt"
    conda:
        "../../envs/samtools.yaml"
    shell:
        "samtools faidx {input} &> {log}"

def get_mummerplot_scaffolds(wildcards):
    filelist = []
    if singularity_enabled():
        for asmname in get_all_accessions():
            if not has_assembly_location(asmname):
                reference = get_reference_id(asmname)
                filelist.append(f"results/{asmname}/2.scaffolding/03.mummer/{asmname}.vs.{reference}.plot.gp")
                filelist.append(f"results/{asmname}/2.scaffolding/03.mummer/{asmname}.vs.{reference}.plot.large.gp")
    return filelist

def get_hic_plots(wildcards):
    filelist = []
    if singularity_enabled():
        for asmname in get_all_accessions():
            if has_hic(asmname) and not has_assembly_location(asmname):
                filelist.append(f"results/{asmname}/2.scaffolding/01.{config['scaffolder']}/contact_map.pdf")
    return filelist

rule scaffold:
    input:
        expand("final_output/{asmname}.full.fa",
            asmname = get_all_accessions()),
        get_mummerplot_scaffolds,
        get_hic_plots,
    output:
        touch("results/scaffolding.done")
