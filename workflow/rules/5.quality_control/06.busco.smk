rule busco_download:
    output:
        temporary(directory("busco_downloads")), #unfortunately cannot be changed as the `--download_path` parameter from busco doesn't work
    log:
        "results/logs/5.quality_control/busco_download.log"
    benchmark:
        "results/benchmarks/5.quality_control/busco_download.txt"
    params:
        odb = config["odb"]
    conda:
        "../../envs/busco.yaml"
    shell:
        "busco --download {params.odb} --download_path {output} &> {log}"

rule compleasm_download:
    output:
        temporary(directory("results/mb_downloads")),
    log:
        "results/logs/5.quality_control/compleasm_download.log"
    benchmark:
        "results/benchmarks/5.quality_control/compleasm_download.txt"
    params:
        lineage = config["odb"].split('_')[0],
    conda:
        "../../envs/busco.yaml"
    shell:
        "compleasm download {params.lineage} -L {output} &> {log}"

rule busco_proteome:
    input:
        proteome = "results/{asmname}/5.quality_control/proteome.pep.fa",
        download = rules.busco_download.output,
    output:
        "results/{asmname}/5.quality_control/06.busco_proteome/{asmname}_proteome/short_summary.specific.{odb}.{asmname}_proteome.txt",
    log:
        "results/logs/5.quality_control/busco_proteome/{odb}/{asmname}.log"
    benchmark:
        "results/benchmarks/5.quality_control/busco_proteome/{odb}/{asmname}.txt"
    threads:
        10
    conda:
        "../../envs/busco.yaml"
    shell:
        "busco -f -i {input.proteome} -l {wildcards.odb} -o {wildcards.asmname}_proteome --out_path $(dirname $(dirname {output})) -m proteome -c {threads} --tar --offline &> {log}"

rule compleasm:
    input:
        genome = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa",
        downloaddir = rules.compleasm_download.output,
    output:
        summary = "results/{asmname}/5.quality_control/06.compleasm/{asmname}/summary.txt",
        hmmer_output = temporary(directory(expand("results/{{asmname}}/5.quality_control/06.compleasm/{{asmname}}/{odb}/hmmer_output", odb=config["odb"]))), #needed because compleasm doesn't compress this
    log:
        "results/logs/5.quality_control/compleasm/{asmname}.log"
    benchmark:
        "results/benchmarks/5.quality_control/compleasm/{asmname}.txt"
    params:
        lineage = config["odb"].split('_')[0],
    threads:
        10
    conda:
        "../../envs/busco.yaml"
    shell:
        "compleasm run -a {input.genome} -o $(dirname {output.summary}) -t {threads} -l {params.lineage} -L {input.downloaddir} &> {log}"

rule compleasm_summary_reformat:
    input:
        "results/{asmname}/5.quality_control/06.compleasm/{asmname}/summary.txt"
    output:
        "results/{asmname}/5.quality_control/06.compleasm/{asmname}/short_summary.specific.{odb}.{asmname}_genome.txt",
    log:
        "results/logs/5.quality_control/compleasm_summary_reformat/{odb}/{asmname}.log"
    benchmark:
        "results/benchmarks/5.quality_control/compleasm_summary_reformat/{odb}/{asmname}.txt"
    shell:
        "awk 'BEGIN{{OFS = \"\\t\";}} /^S/{{busco[\"Complete and single-copy BUSCOs\"]=$2;}} /^D/{{busco[\"Complete and duplicated BUSCOs\"]=$2;}} /^[FI]/{{busco[\"Fragmented BUSCOs\"]+=$2;}} /^M/{{busco[\"Missing BUSCOs\"]=$2;}} END{{for (i in busco){{print i,busco[i];}}}}' {input} > {output} 2> {log}"

rule busco_plot:
    input:
        lambda wildcards: expand("results/{asmname}/5.quality_control/06.compleasm/{asmname}/short_summary.specific.{odb}.{asmname}_genome.txt", asmname=config["set"][wildcards.asmset], odb=config["odb"]),
        lambda wildcards: expand("results/{asmname}/5.quality_control/06.busco_proteome/{asmname}_proteome/short_summary.specific.{odb}.{asmname}_proteome.txt", asmname=config["set"][wildcards.asmset], odb=config["odb"]),
    output:
        report("results/{asmset}/5.quality_control/06.busco_plot/busco_figure.png", category="Gene completeness", labels={"type": "busco", "set": "{asmset}"}),
    log:
        "results/logs/5.quality_control/busco_plot/{asmset}.log"
    benchmark:
        "results/benchmarks/5.quality_control/busco_plot/{asmset}.txt"
    conda:
        "../../envs/busco.yaml"
    shell:
        """
        (
        cp {input} $(dirname {output})/
        generate_plot.py -wd $(dirname {output})
        ) &> {log}
        """