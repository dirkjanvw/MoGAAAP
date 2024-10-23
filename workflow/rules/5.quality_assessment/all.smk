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
    k = config["k_qa"]
    all_output = []

    # HiFi
    for asmname in get_all_accessions():
        all_output.append(f"results/{asmname}/5.quality_assessment/01.merqury/{k}/hifi/{asmname}_vs_hifi.{asmname}.qv.html")  #per sequence qv
        all_output.append(f"results/{asmname}/5.quality_assessment/01.merqury/{k}/hifi/{asmname}_vs_hifi.{asmname}.spectra-cn.fl.png")  #spectra-cn

    # ONT (optional)
    for asmname in get_all_accessions():
        if not SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(asmname)]["ont"].isnull().values.item():
            all_output.append(f"results/{asmname}/5.quality_assessment/01.merqury/{k}/ont/{asmname}_vs_ont.{asmname}.qv.html")  #per sequence qv
            all_output.append(f"results/{asmname}/5.quality_assessment/01.merqury/{k}/ont/{asmname}_vs_ont.{asmname}.spectra-cn.fl.png")  #spectra-cn

    # Illumina (optional)
    for asmname in get_all_accessions():
        if not SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(asmname)]["illumina_1"].isnull().values.item():
            all_output.append(f"results/{asmname}/5.quality_assessment/01.merqury/{k}/illumina/{asmname}_vs_illumina.{asmname}.qv.html")  #per sequence qv
            all_output.append(f"results/{asmname}/5.quality_assessment/01.merqury/{k}/illumina/{asmname}_vs_illumina.{asmname}.spectra-cn.fl.png")  #spectra-cn

    return all_output

def get_multiqc_output(wildcards):
    all_output = []
    for asmname in get_all_accessions():
        if not SAMPLES[SAMPLES["accessionId"] == get_clean_accession_id(asmname)]["illumina_1"].isnull().values.item():
            all_output.append(f"results/{asmname}/5.quality_assessment/04.multiqc/multiqc_report.html")  #mapping
    return all_output

def get_pantools_output(wildcards):
    all_output = []
    for asmset in config["set"]:
        if (len(get_all_accessions_from_asmset(asmset)) < 2):
            continue
        all_output.append(f"results/{asmset}/5.quality_assessment/05.pantools/panproteome_groups_DB/pantools_homology_groups.txt")  #pantools homology grouping
        if len(config["set"][asmset]) >= 3:
            all_output.append(f"results/{asmset}/5.quality_assessment/05.pantools/panproteome_groups_DB/pangenome_size/gene/core_dispensable_growth.png")  #pantools pangenome growth
            if len(config["set"][asmset]) <= 10:
                    all_output.append(f"results/{asmset}/5.quality_assessment/05.pantools/panproteome_groups_DB/gene_classification/upset/output/genomes.pdf")  #pantools gene classification
    return all_output

def get_busco_output(wildcards):
    all_output = []
    for asmset in config["set"]:
        all_output.append(f"results/{asmset}/5.quality_assessment/06.busco_plot/busco_figure.png")  #busco (proteome) and compleasm (genome)
    return all_output

def get_omark_output(wildcards):
    all_output = []
    for asmset in config["set"]:
        all_output.append(f"results/{asmset}/5.quality_assessment/07.omark_plot.png")  #omark
    return all_output

def get_kmerdb_output(wildcards):
    all_output = []
    k = config["k_qa"]
    for asmset in config["set"]:
        if (len(get_all_accessions_from_asmset(asmset)) < 2):
            continue
        if len(config["set"][asmset]) >= 2:
            all_output.append(f"results/{asmset}/5.quality_assessment/08.kmer-db/{k}/{asmset}.csv.mash.pdf"),  #kmer distances
    return all_output

def get_mash_output(wildcards):
    all_output = []
    for asmset in config["set"]:
        if (len(get_all_accessions_from_asmset(asmset)) < 2):
            continue
        if len(config["set"][asmset]) >= 2:
            all_output.append(f"results/{asmset}/5.quality_assessment/09.mash/{asmset}.pdf"),  #mash distances
    return all_output

def get_ntsynt_output(wildcards):
    all_output = []
    mink = 24
    minw = 1000
    for asmset in config["set"]:
        if (len(get_all_accessions_from_asmset(asmset)) < 2):
            continue
        if len(config["set"][asmset]) >= 2:  #synteny only makes sense for multiple genomes
            all_output.append(f"results/{asmset}/5.quality_assessment/10.ntsynt/{asmset}.k{mink}.w{minw}.png")
    return all_output

def get_sans_output(wildcards):
    all_output = []
    k = config["k_qa"]
    bootstrap = 1000
    for asmset in config["set"]:
        if (len(get_all_accessions_from_asmset(asmset)) < 2):
            continue
        if len(config["set"][asmset]) >= 4:  #phylogenetic networks only makes sense for 4+ genomes
            all_output.append(f"results/{asmset}/5.quality_assessment/11.sans/{k}/{asmset}_b{bootstrap}.nexus")
    return all_output

def get_pangrowth_output(wildcards):
    all_output = []
    k = config["k_qa"]
    for asmset in config["set"]:
        if (len(get_all_accessions_from_asmset(asmset)) < 2):
            continue
        if len(config["set"][asmset]) >= 3:
            all_output.append(f"results/{asmset}/5.quality_assessment/12.pangrowth/{k}/hist.pdf"),  #pangrowth hist
            all_output.append(f"results/{asmset}/5.quality_assessment/12.pangrowth/{k}/growth.pdf"),  #pangrowth growth
            all_output.append(f"results/{asmset}/5.quality_assessment/12.pangrowth/{k}/core.pdf"),  #pangrowth core
    return all_output

def get_statistics_output(wildcards):
    all_output = []
    for asmset in config["set"]:
        all_output.append(f"results/{asmset}/5.quality_assessment/13.statistics/{asmset}.html")
    return all_output

rule qa:
    input:
        # individual outputs
        get_merqury_output,  #merqury
        expand("results/{asmname}/5.quality_assessment/02.kraken2/{asmname}.kraken2.krona.html", asmname=get_all_accessions()),  #kraken2
        expand("results/{asmname}/5.quality_assessment/03.fcs/{asmname}.fcs_gx_report.html", asmname=get_all_accessions()),  #fcs-gx
        expand("results/{asmname}/5.quality_assessment/03.fcs/{asmname}/fcs_adaptor_report.html", asmname=get_all_accessions()),  #fcs-adaptor
        get_multiqc_output,  #mapping

        # grouped outputs
        get_pantools_output,  #pantools
        get_busco_output,  #busco
        get_omark_output,  #omark
        get_kmerdb_output,  #kmer-db
        get_mash_output,  #mash
        get_ntsynt_output, #ntsynt
        get_sans_output,  #sans nexus file (genome only with 1000 bootstrap)
        get_pangrowth_output,  #pangrowth
        get_statistics_output, #statistics
    output:
        touch("results/quality_assessment.done")
