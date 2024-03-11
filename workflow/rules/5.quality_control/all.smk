include: "01.merqury.smk"
include: "02.kraken.smk"
include: "03.fcs.smk"
include: "04.pantools.smk"
include: "05.busco.smk"
include: "06.omark.smk"
include: "07.kmer-db.smk"
include: "08.mash.smk"
include: "09.ntsynt.smk"
include: "10.sans.smk"
include: "11.pangrowth.smk"
include: "12.collinearity.smk"
include: "13.mapping.smk"

def get_merqury_output(wildcards):
    k = config["k"]
    all_output = []

    # HiFi
    for asmname in config["reads"]["hifi"]:
        all_output.append(f"results/{asmname}/5.quality_control/01.merqury/{k}/hifi/{asmname}_vs_hifi.qv")  #merqury
        all_output.append(f"results/{asmname}/5.quality_control/01.merqury/{k}/hifi/{asmname}_vs_hifi.{asmname}.spectra-cn.fl.png")  #merqury

    # ONT
    for asmname in config["reads"]["ont"]:
        all_output.append(f"results/{asmname}/5.quality_control/01.merqury/{k}/ont/{asmname}_vs_ont.qv")  #merqury
        all_output.append(f"results/{asmname}/5.quality_control/01.merqury/{k}/ont/{asmname}_vs_ont.{asmname}.spectra-cn.fl.png")  #merqury

    # Illumina
    for asmname in config["reads"]["illumina"]:
        for sample in config["reads"]["illumina"]:
            all_output.append(f"results/{asmname}/5.quality_control/01.merqury/{k}/{sample}/{asmname}_vs_{sample}.qv")  #merqury
            all_output.append(f"results/{asmname}/5.quality_control/01.merqury/{k}/{sample}/{asmname}_vs_{sample}.{asmname}.spectra-cn.fl.png")  #merqury

    return all_output

def get_pantools_output(wildcards):
    all_output = []
    for asmset in config["set"]:
        all_output.append(f"results/{asmset}/5.quality_control/04.pantools/panproteome_groups_DB/pantools_homology_groups.txt")  #pantools homology grouping
        if len(get_all_accessions()) >= 3:
            all_output.append(f"results/{asmset}/5.quality_control/04.pantools/panproteome_groups_DB/pangenome_size/gene/core_dispensable_growth.png")  #pantools pangenome growth
            if len(get_all_accessions()) <= 10:
                    all_output.append(f"results/{asmset}/5.quality_control/04.pantools/panproteome_groups_DB/gene_classification/upset/output/genomes.pdf")  #pantools gene classification
    return all_output

rule qc:
    input:
        # individual outputs
        get_merqury_output,  #merqury
        expand("results/{asmname}/5.quality_control/02.kraken2/{asmname}.kraken2.krona.html", asmname=get_all_accessions()),  #kraken2
        expand("results/{asmname}/5.quality_control/03.fcs/{asmname}.fcs_gx_report.txt", asmname=get_all_accessions()),  #fcs-gx
        expand("results/{asmname}/5.quality_control/03.fcs/{asmname}.fcs_adaptor_report.txt", asmname=get_all_accessions()),  #fcs-adaptor

        # grouped outputs
        get_pantools_output,  #pantools
        expand("results/{asmset}/5.quality_control/05.busco_plot/busco_figure.png", asmset=config["set"]),  #busco (proteome) and compleasm (genome)
        expand("results/{asmset}/5.quality_control/06.omark_plot.png", asmset=config["set"]),  #omark
        expand("results/{asmset}/5.quality_control/07.kmer-db/k{k}.csv.mash.pdf", asmset=config["set"], k=config["k"]),  #kmer distances
        expand("results/{asmset}/5.quality_control/08.mash/{asmset}.pdf", asmset=config["set"]), #mash distances
        expand("results/{asmset}/5.quality_control/09.ntsynt/{asmset}.k{mink}.w{minw}.png", asmset=config["set"], mink=24, minw=1000), #ntsynt
        expand("results/{asmset}/5.quality_control/10.sans/{k}/{asmset}_b{bootstrap}.nexus", k=config["k"], asmset=config["set"], bootstrap=1000),  #sans nexus file (genome only with 1000 bootstrap)
        expand("results/{asmset}/5.quality_control/11.pangrowth/{k}/{figure}.pdf", asmset=config["set"], k=config["k"], figure=["hist", "growth", "core"]),  #pangrowth
        expand("results/{asmset}/5.quality_control/12.collinearity/karyotype.pdf", asmset=config["set"]),  #jcvi collinearity
        expand("results/{asmset}/5.quality_control/13.multiqc/multiqc_report.html", asmset=config["set"]),  #mapping
    output:
        touch("results/quality_control.done")
