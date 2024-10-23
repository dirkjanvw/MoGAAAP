rule omamer_search:
    input:
        "results/{asmname}/5.quality_assessment/proteome.pep.fa",
    output:
        "results/{asmname}/5.quality_assessment/07.omamer/{asmname}.omamer",
    log:
        "results/logs/5.quality_assessment/omamer/{asmname}.log"
    benchmark:
        "results/benchmarks/5.quality_assessment/omamer/{asmname}.txt"
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
        omamer = "results/{asmname}/5.quality_assessment/07.omamer/{asmname}.omamer",
        splice = "results/{asmname}/5.quality_assessment/proteome.splice",
    output:
        omamer_omq = "results/{asmname}/5.quality_assessment/07.omark/{asmname}.omq",
        omamer_pdf = "results/{asmname}/5.quality_assessment/07.omark/{asmname}.pdf",
        omamer_png = "results/{asmname}/5.quality_assessment/07.omark/{asmname}.png",
        omamer_sum = "results/{asmname}/5.quality_assessment/07.omark/{asmname}.sum",
        omamer_tax = "results/{asmname}/5.quality_assessment/07.omark/{asmname}.tax",
        omamer_ump = "results/{asmname}/5.quality_assessment/07.omark/{asmname}.ump",
        omamer_sumtxt = "results/{asmname}/5.quality_assessment/07.omark/{asmname}_detailed_summary.txt",
        omamer_isotxt = "results/{asmname}/5.quality_assessment/07.omark/{asmname}_selected_isoforms.txt",
    log:
        "results/logs/5.quality_assessment/omark/{asmname}.log"
    benchmark:
        "results/benchmarks/5.quality_assessment/omark/{asmname}.txt"
    params:
        database = config["OMAdb"]
    conda:
        "../../envs/oma.yaml"
    shell:
        "omark -f {input.omamer} -i {input.splice} -d {params.database} -o $(dirname {output.omamer_sum}) -v &> {log}"

rule omark_plot:
    input:
        lambda wildcards: expand("results/{asmname}/5.quality_assessment/07.omark/{asmname}.sum", asmname=get_all_accessions_from_asmset(wildcards.asmset, 1)),
    output:
        tmpdir = temporary(directory("results/{asmset}/5.quality_assessment/07.omark_plot")),
        png = report("results/{asmset}/5.quality_assessment/07.omark_plot.png",
            category="Gene completeness",
            caption="../../report/omark.rst",
            labels={"type": "omark", "set": "{asmset}"}),
    log:
        "results/logs/5.quality_assessment/omark_plot/{asmset}.log"
    benchmark:
        "results/benchmarks/5.quality_assessment/omark_plot/{asmset}.txt"
    conda:
        "../../envs/oma.yaml"
    shell:
        """
        (
        mkdir {output.tmpdir}
        cp {input} {output.tmpdir}
        workflow/scripts/OMArk.a3c75ad/plot_all_results.py -i $(dirname {output.tmpdir}) -o {output.png}
        ) &> {log}
        """
