rule busco_download:
    output:
        temporary(directory("busco_downloads")), #unfortunately cannot be changed as the `--download_path` parameter from busco doesn't work
    log:
        "results/logs/3.quality_assessment/busco_download.log"
    benchmark:
        "results/benchmarks/3.quality_assessment/busco_download.txt"
    params:
        odb = config["odb"]
    conda:
        "../../envs/busco.yaml"
    shell:
        "busco --download {params.odb} --download_path {output} &> {log}"

rule busco_proteome:
    input:
        proteome = "results/{asmname}/3.quality_assessment/proteome.pep.fa",
        download = ancient(rules.busco_download.output),
    output:
        "results/{asmname}/3.quality_assessment/06.busco_proteome/{asmbase}_proteome/short_summary.specific.{odb}.{asmbase}_proteome.txt",
    log:
        "results/logs/3.quality_assessment/busco_proteome/{odb}/{asmname}/{asmbase}.log"
    benchmark:
        "results/benchmarks/3.quality_assessment/busco_proteome/{odb}/{asmname}/{asmbase}.txt"
    threads:
        10
    conda:
        "../../envs/busco.yaml"
    shell:
        "busco -f -i {input.proteome} -l {wildcards.odb} -o $(basename $(dirname {output})) --out_path $(dirname $(dirname {output})) -m proteome -c {threads} --tar --offline &> {log}"

rule busco_genome:
    input:
        genome = "final_output/{asmname}.full.fa",
        download = ancient(rules.busco_download.output),
    output:
        "results/{asmname}/3.quality_assessment/06.busco_genome/{asmbase}_genome/short_summary.specific.{odb}.{asmbase}_genome.txt",
    log:
        "results/logs/3.quality_assessment/busco_genome/{odb}/{asmname}/{asmbase}.log"
    benchmark:
        "results/benchmarks/3.quality_assessment/busco_genome/{odb}/{asmname}/{asmbase}.txt"
    threads:
        10
    conda:
        "../../envs/busco.yaml"
    shell:
        "busco -f -i {input.genome} -l {wildcards.odb} -o $(basename $(dirname {output})) --out_path $(dirname $(dirname {output})) -m genome -c {threads} --tar --offline &> {log}"

def get_busco_plot_input(wildcards):
    odb = config["odb"]
    all_output = []
    for asmname in get_all_accessions_from_asmset(wildcards.asmset, 1):
        if get_haplotype_information(asmname) > 1:
            asmbase = get_clean_accession_id(asmname)
            haplotype = get_haplotype_accession_id(asmname)
            all_output.append(f"results/{asmbase}.{haplotype}/3.quality_assessment/06.busco_genome/{asmbase}_{haplotype}_genome/short_summary.specific.{odb}.{asmbase}_{haplotype}_genome.txt")
            all_output.append(f"results/{asmbase}.{haplotype}/3.quality_assessment/06.busco_proteome/{asmbase}_{haplotype}_proteome/short_summary.specific.{odb}.{asmbase}_{haplotype}_proteome.txt")
        else:
            all_output.append(f"results/{asmname}/3.quality_assessment/06.busco_genome/{asmname}_genome/short_summary.specific.{odb}.{asmname}_genome.txt")
            all_output.append(f"results/{asmname}/3.quality_assessment/06.busco_proteome/{asmname}_proteome/short_summary.specific.{odb}.{asmname}_proteome.txt")
    return all_output

rule busco_plot:
    input:
        get_busco_plot_input,
    output:
        report("results/{asmset}/3.quality_assessment/06.busco_plot/busco_figure.png",
            category="Quality assessment",
            subcategory="Gene completeness",
            caption="../../report/busco.rst",
            labels={"type": "busco", "set": "{asmset}"}),
    log:
        "results/logs/3.quality_assessment/busco_plot/{asmset}.log"
    benchmark:
        "results/benchmarks/3.quality_assessment/busco_plot/{asmset}.txt"
    conda:
        "../../envs/busco.yaml"
    shell:
        """
        (
        cp {input} $(dirname {output})/
        generate_plot.py -wd $(dirname {output})
        ) &> {log}
        """
