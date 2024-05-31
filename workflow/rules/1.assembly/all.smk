include: "01.hifiasm.smk"
include: "02.contigs.smk"
include: "03.mummer.smk"

def get_mummerplot_contigs(wildcards):
    filelist = []
    minlen = config["min_contig_len"]
    for asmname in get_all_accessions():
        reference = get_reference_id(asmname)
        filelist.append(f"results/{asmname}/1.assembly/03.mummer/{asmname}.min{minlen}.vs.{reference}.plot.gp")
    return filelist

rule assemble:
    input:
        expand("results/{asmname}/1.assembly/02.contigs/{asmname}.min{minlen}.sorted.renamed.fa",
            asmname=get_all_accessions(),
            minlen=config["min_contig_len"]
        ),
        get_mummerplot_contigs,
    output:
        touch("results/assembly.done")
