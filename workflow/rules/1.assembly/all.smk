include: "01.assembly.smk"
include: "02.contigs.smk"
include: "03.mummer.smk"

rule copy_contigs:
    input:
        expand("results/{{asmname}}/1.assembly/02.contigs/{{asmname}}.min{minlen}.sorted.renamed.fa", minlen=config["min_contig_len"]),
    output:
        protected("final_output/{asmname}.contigs.fa"),
    log:
        "results/logs/2.scaffolding/copy_contigs/{asmname}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/copy_contigs/{asmname}.txt"
    shell:
        "cp $(realpath {input}) {output} &> {log}"

def get_mummerplot_contigs(wildcards):
    filelist = []
    minlen = config["min_contig_len"]
    for asmname in get_all_accessions():
        reference = get_reference_id(asmname)
        filelist.append(f"results/{asmname}/1.assembly/03.mummer/{asmname}.min{minlen}.vs.{reference}.plot.gp")
    return filelist

rule assemble:
    input:
        expand("final_output/{asmname}.contigs.fa",
            asmname=get_all_accessions()
        ),
        get_mummerplot_contigs,
    output:
        touch("results/assembly.done")
