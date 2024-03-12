rule omamer_search:
    input:
        "results/{asmname}/5.quality_control/proteome.pep.fa",
    output:
        "results/{asmname}/5.quality_control/07.omamer/{asmname}.omamer",
    log:
        "results/logs/5.quality_control/omamer/{asmname}.log"
    benchmark:
        "results/benchmarks/5.quality_control/omamer/{asmname}.txt"
    params:
        database = config["OMAdb"]
    threads:
        10
    conda:
        "../../envs/oma.yaml"
    shell:
        "omamer search -d {params.database} -q {input} -o {output} -t {threads} &> {log}"

rule omark:
    input:
        omamer = "results/{asmname}/5.quality_control/07.omamer/{asmname}.omamer",
        splice = "results/{asmname}/5.quality_control/transcriptome.splice",
    output:
        omamer_omq = "results/{asmname}/5.quality_control/07.omark/{asmname}.omq",
        omamer_pdf = "results/{asmname}/5.quality_control/07.omark/{asmname}.pdf",
        omamer_png = "results/{asmname}/5.quality_control/07.omark/{asmname}.png",
        omamer_sum = "results/{asmname}/5.quality_control/07.omark/{asmname}.sum",
        omamer_tax = "results/{asmname}/5.quality_control/07.omark/{asmname}.tax",
        omamer_ump = "results/{asmname}/5.quality_control/07.omark/{asmname}.ump",
        omamer_sumtxt = "results/{asmname}/5.quality_control/07.omark/{asmname}_detailed_summary.txt",
        omamer_isotxt = "results/{asmname}/5.quality_control/07.omark/{asmname}_selected_isoforms.txt",
    log:
        "results/logs/5.quality_control/omark/{asmname}.log"
    benchmark:
        "results/benchmarks/5.quality_control/omark/{asmname}.txt"
    params:
        database = config["OMAdb"]
    conda:
        "../../envs/oma.yaml"
    shell:
        "omark -f {input.omamer} -i {input.splice} -d {params.database} -o $(dirname {output.omamer_sum}) -v &> {log}"

rule omark_plot:
    input:
        lambda wildcards: expand("results/{asmname}/5.quality_control/07.omark/{asmname}.sum", asmname=config["set"][wildcards.asmset]),
    output:
        report("results/{asmset}/5.quality_control/07.omark_plot.png", category="Gene completeness", labels={"type": "omark", "set": "{asmset}"}),
    log:
        "results/logs/5.quality_control/omark_plot/{asmset}.log"
    benchmark:
        "results/benchmarks/5.quality_control/omark_plot/{asmset}.txt"
    conda:
        "../../envs/oma.yaml"
    shell:
        "workflow/scripts/OMArk.a3c75ad/plot_all_results.py -i $(dirname $(dirname {input}) | uniq) -o {output} &> {log}"