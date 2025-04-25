include: "01.assembly.smk"
include: "02.contigs.smk"
include: "03.mummer.smk"
include: "04.scaffolding.smk"
include: "05.renaming.smk"
include: "06.mummer.smk"

rule copy_contigs:
    input:
        expand("results/{{asmname}}/1.assembly/02.contigs/{{asmname}}.min{minlen}.sorted.renamed.fa", minlen=config["min_contig_len"]),
    output:
        "final_output/{asmname}.contigs.fa",
    log:
        "results/logs/1.assembly/copy_contigs/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/copy_contigs/{asmname}.txt"
    shell:
        "cp $(realpath {input}) {output} &> {log}"

def get_mummerplot_contigs(wildcards):
    filelist = []
    if singularity_enabled():
        minlen = config["min_contig_len"]
        for asmname in get_all_accessions():
            if not has_assembly_location(asmname):
                reference = get_reference_id(asmname)
                filelist.append(f"results/{asmname}/1.assembly/03.mummer/{asmname}.min{minlen}.vs.{reference}.plot.gp")
                filelist.append(f"results/{asmname}/1.assembly/03.mummer/{asmname}.min{minlen}.vs.{reference}.plot.large.gp")
    return filelist

rule copy_assembly:
    input:
        lambda wildcards:
            "results/{asmname}/1.assembly/05.renaming/{asmname}.fa"
            if not has_assembly_location(wildcards.asmname)
            else get_assembly_location(wildcards.asmname),
    output:
        "final_output/{asmname}.full.fa",
    log:
        "results/logs/1.assembly/copy_assembly/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/copy_assembly/{asmname}.txt"
    conda:
        "../../envs/seqkit.yaml"
    shell:
        "seqkit seq -w 80 $(realpath {input}) -o {output} &> {log}"

rule index_scaffolds:
    input:
        "final_output/{asmname}.full.fa"
    output:
        "final_output/{asmname}.full.fa.fai"
    log:
        "results/logs/1.assembly/index_scaffolds/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/index_scaffolds/{asmname}.txt"
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
                filelist.append(f"results/{asmname}/1.assembly/06.mummer/{asmname}.vs.{reference}.plot.gp")
                filelist.append(f"results/{asmname}/1.assembly/06.mummer/{asmname}.vs.{reference}.plot.large.gp")
    return filelist

def get_hic_plots(wildcards):
    filelist = []
    if singularity_enabled():
        for asmname in get_all_accessions():
            if has_hic(asmname) and not has_assembly_location(asmname): #TODO: enable for all assemblies
                filelist.append(f"results/{asmname}/1.assembly/04.{config['scaffolder']}/contact_map.pdf")
    return filelist

rule contig:
    input:
        expand("final_output/{asmname}.contigs.fa",
            asmname=[asmname for asmname in get_all_accessions() if not has_assembly_location(asmname)]
        ),
        get_mummerplot_contigs,
    output:
        touch("results/contig.done")

rule assemble:
    input:
        "results/contig.done",
        expand("final_output/{asmname}.full.fa",
            asmname=[asmname for asmname in get_all_accessions() if not has_assembly_location(asmname)]
        ),
        get_mummerplot_scaffolds,
        get_hic_plots,
    output:
        touch("results/assembly.done")
