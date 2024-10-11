rule busco_download:
    output:
        temporary(directory("busco_downloads")), #unfortunately cannot be changed as the `--download_path` parameter from busco doesn't work
    log:
        "results/logs/5.quality_assessment/busco_download.log"
    benchmark:
        "results/benchmarks/5.quality_assessment/busco_download.txt"
    params:
        odb = config["odb"]
    conda:
        "../../envs/busco.yaml"
    shell:
        "busco --download {params.odb} --download_path {output} &> {log}"

rule busco_proteome:
    input:
        proteome = "results/{asmname}/5.quality_assessment/proteome.pep.fa",
        download = rules.busco_download.output,
    output:
        "results/{asmname}/5.quality_assessment/06.busco_proteome/{asmname}_proteome/short_summary.specific.{odb}.{asmname}_proteome.txt",
    log:
        "results/logs/5.quality_assessment/busco_proteome/{odb}/{asmname}.log"
    benchmark:
        "results/benchmarks/5.quality_assessment/busco_proteome/{odb}/{asmname}.txt"
    threads:
        10
    conda:
        "../../envs/busco.yaml"
    shell:
        "busco -f -i {input.proteome} -l {wildcards.odb} -o {wildcards.asmname}_proteome --out_path $(dirname $(dirname {output})) -m proteome -c {threads} --tar --offline &> {log}"

rule busco_genome:
    input:
        genome = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa",
        download = rules.busco_download.output,
    output:
        "results/{asmname}/5.quality_assessment/06.busco_genome/{asmname}_genome/short_summary.specific.{odb}.{asmname}_genome.txt",
    log:
        "results/logs/5.quality_assessment/busco_genome/{odb}/{asmname}.log"
    benchmark:
        "results/benchmarks/5.quality_assessment/busco_genome/{odb}/{asmname}.txt"
    threads:
        10
    conda:
        "../../envs/busco.yaml"
    shell:
        "busco -f -i {input.genome} -l {wildcards.odb} -o {wildcards.asmname}_genome --out_path $(dirname $(dirname {output})) -m genome -c {threads} --tar --offline &> {log}"

rule busco_plot:
    input:
        lambda wildcards: expand("results/{asmname}/5.quality_assessment/06.busco_genome/{asmname}_genome/short_summary.specific.{odb}.{asmname}_genome.txt", asmname=get_all_accessions_from_asmset(wildcards.asmset), odb=config["odb"]),
        lambda wildcards: expand("results/{asmname}/5.quality_assessment/06.busco_proteome/{asmname}_proteome/short_summary.specific.{odb}.{asmname}_proteome.txt", asmname=get_all_accessions_from_asmset(wildcards.asmset), odb=config["odb"]),
    output:
        report("results/{asmset}/5.quality_assessment/06.busco_plot/busco_figure.png",
            category="Gene completeness",
            caption="../../report/busco.rst",
            labels={"type": "busco", "set": "{asmset}"}),
    log:
        "results/logs/5.quality_assessment/busco_plot/{asmset}.log"
    benchmark:
        "results/benchmarks/5.quality_assessment/busco_plot/{asmset}.txt"
    conda:
        "../../envs/busco.yaml"
    shell:
        """
        (
        cp {input} $(dirname {output})/
        generate_plot.py -wd $(dirname {output})
        ) &> {log}
        """