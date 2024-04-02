include: "00.proteome.smk"
include: "01.merqury.smk"
include: "02.kraken.smk"
include: "03.fcs.smk"
include: "04.mapping.smk"
include: "05.pantools.smk"
include: "06.busco.smk"
include: "07.omark.smk"
include: "08.kmer-db.smk"
include: "09.mash.smk"
include: "10.ntsynt.smk"
include: "11.sans.smk"
include: "12.pangrowth.smk"
include: "13.statistics.smk"

def get_merqury_output(wildcards):
    k = config["k_qc"]
    all_output = []

    # HiFi
    for asmname in config["reads"]["hifi"]:
        for sample in config["reads"]["hifi"][asmname]:
            all_output.append(f"results/{asmname}/5.quality_control/01.merqury/{k}/{sample}/{asmname}_vs_{sample}.{asmname}.qv.pdf")  #per sequence qv
            all_output.append(f"results/{asmname}/5.quality_control/01.merqury/{k}/{sample}/{asmname}_vs_{sample}.{asmname}.spectra-cn.fl.png")  #spectra-cn

    # ONT (optional)
    if "ont" in config["reads"]:
        for asmname in config["reads"]["ont"]:
            for sample in config["reads"]["ont"][asmname]:
                all_output.append(f"results/{asmname}/5.quality_control/01.merqury/{k}/{sample}/{asmname}_vs_{sample}.{asmname}.qv.pdf")  #per sequence qv
                all_output.append(f"results/{asmname}/5.quality_control/01.merqury/{k}/{sample}/{asmname}_vs_{sample}.{asmname}.spectra-cn.fl.png")  #spectra-cn

    # Illumina (optional)
    if "illumina" in config["reads"]:
        for asmname in config["reads"]["illumina"]:
            for sample in config["reads"]["illumina"][asmname]:
                all_output.append(f"results/{asmname}/5.quality_control/01.merqury/{k}/{sample}/{asmname}_vs_{sample}.{asmname}.qv.pdf")  #per sequence qv
                all_output.append(f"results/{asmname}/5.quality_control/01.merqury/{k}/{sample}/{asmname}_vs_{sample}.{asmname}.spectra-cn.fl.png")  #spectra-cn

    return all_output

def get_multiqc_output(wildcards):
    all_output = []
    if "illumina" in config["reads"]:
        for asmname in get_all_accessions():
            if asmname in config["reads"]["illumina"]:
                all_output.append(f"results/{asmname}/5.quality_control/04.multiqc/multiqc_report.html")  #mapping
    return all_output

def get_pantools_output(wildcards):
    all_output = []
    for asmset in config["set"]:
        all_output.append(f"results/{asmset}/5.quality_control/05.pantools/panproteome_groups_DB/pantools_homology_groups.txt")  #pantools homology grouping
        if len(config["set"][asmset]) >= 3:
            all_output.append(f"results/{asmset}/5.quality_control/05.pantools/panproteome_groups_DB/pangenome_size/gene/core_dispensable_growth.png")  #pantools pangenome growth
            if len(config["set"][asmset]) >= 10:
                    all_output.append(f"results/{asmset}/5.quality_control/05.pantools/panproteome_groups_DB/gene_classification/upset/output/genomes.pdf")  #pantools gene classification
    return all_output

def get_pangrowth_output(wildcards):
    all_output = []
    k = config["k_qc"]
    for asmset in config["set"]:
        if len(config["set"][asmset]) >= 3:
            all_output.append(f"results/{asmset}/5.quality_control/12.pangrowth/{k}/hist.pdf"),  #pangrowth hist
            all_output.append(f"results/{asmset}/5.quality_control/12.pangrowth/{k}/growth.pdf"),  #pangrowth growth
            all_output.append(f"results/{asmset}/5.quality_control/12.pangrowth/{k}/core.pdf"),  #pangrowth core
    return all_output

rule qc:
    input:
        # individual outputs
        get_merqury_output,  #merqury
        expand("results/{asmname}/5.quality_control/02.kraken2/{asmname}.kraken2.krona.html", asmname=get_all_accessions()),  #kraken2
        expand("results/{asmname}/5.quality_control/03.fcs/{asmname}.fcs_gx_report.pdf", asmname=get_all_accessions()),  #fcs-gx
        expand("results/{asmname}/5.quality_control/03.fcs/{asmname}/fcs_adaptor_report.pdf", asmname=get_all_accessions()),  #fcs-adaptor
        get_multiqc_output,  #mapping

        # grouped outputs
        get_pantools_output,  #pantools
        expand("results/{asmset}/5.quality_control/06.busco_plot/busco_figure.png", asmset=config["set"]),  #busco (proteome) and compleasm (genome)
        expand("results/{asmset}/5.quality_control/07.omark_plot.png", asmset=config["set"]),  #omark
        expand("results/{asmset}/5.quality_control/08.kmer-db/{k}/{asmset}.csv.mash.pdf", asmset=config["set"], k=config["k_qc"]),  #kmer distances
        expand("results/{asmset}/5.quality_control/09.mash/{asmset}.pdf", asmset=config["set"]), #mash distances
        expand("results/{asmset}/5.quality_control/10.ntsynt/{asmset}.k{mink}.w{minw}.png", asmset=config["set"], mink=24, minw=1000), #ntsynt
        expand("results/{asmset}/5.quality_control/11.sans/{k}/{asmset}_b{bootstrap}.nexus", k=config["k_qc"], asmset=config["set"], bootstrap=1000),  #sans nexus file (genome only with 1000 bootstrap)
        get_pangrowth_output,  #pangrowth
        expand("results/{asmset}/5.quality_control/13.statistics/{asmset}.tsv", asmset=config["set"]),  #statistics
    output:
        touch("results/quality_control.done")
