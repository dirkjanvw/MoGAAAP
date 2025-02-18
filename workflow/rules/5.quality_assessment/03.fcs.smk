rule fcs_gx:
    input:
        "final_output/{asmname}.full.fa"
    output:
        report = "results/{asmname}/5.quality_assessment/03.fcs/{asmname}.fcs_gx_report.txt",
        taxonomy = "results/{asmname}/5.quality_assessment/03.fcs/{asmname}.taxonomy.rpt",
    log:
        "results/logs/5.quality_assessment/fcs-gx/{asmname}.log"
    benchmark:
        "results/benchmarks/5.quality_assessment/fcs-gx/{asmname}.txt"
    params:
        taxid = lambda wildcards: get_taxid(wildcards),
        gxdb = config["gxdb"],
    resources:
        gbmem = 500
    container:
        "docker://ncbi/fcs-gx:0.5.4"
    retries:
        3
    shell:
        "run_gx --fasta {input} --tax-id {params.taxid} --gx-db {params.gxdb} --out-basename {wildcards.asmname} --out-dir $(dirname {output.taxonomy}) &> {log}"

rule visualise_fcs_gx:
    input:
        "results/{asmname}/5.quality_assessment/03.fcs/{asmname}.fcs_gx_report.txt"
    output:
        tsv = "results/{asmname}/5.quality_assessment/03.fcs/{asmname}.fcs_gx_report.tsv",
        html = report("results/{asmname}/5.quality_assessment/03.fcs/{asmname}.fcs_gx_report.html",
            category="Contamination",
            caption="../../report/fcs-gx.rst",
            labels={"type": "fcs-gx", "assembly": "{asmname}"}),
    log:
        "results/logs/5.quality_assessment/visualise_fcs_gx/{asmname}.log"
    benchmark:
        "results/benchmarks/5.quality_assessment/visualise_fcs_gx/{asmname}.txt"
    conda:
        "../../envs/csvtotable.yaml"
    shell:
        """
        (
        tail -n+2 {input} > {output.tsv}
        csvtotable -d $'\\t' {output.tsv} {output.html}
        ) &> {log}
        """

rule fcs_adaptor:
    input:
        "final_output/{asmname}.full.fa"
    output:
        cleaned = "results/{asmname}/5.quality_assessment/03.fcs/{asmname}/cleaned_sequences/{asmname}.full.fa",
        calls = "results/{asmname}/5.quality_assessment/03.fcs/{asmname}/combined.calls.jsonl",
        fcslog = "results/{asmname}/5.quality_assessment/03.fcs/{asmname}/fcs.log",
        fcsadaptorlog = "results/{asmname}/5.quality_assessment/03.fcs/{asmname}/fcs_adaptor.log",
        fcsadaptorrpt = "results/{asmname}/5.quality_assessment/03.fcs/{asmname}/fcs_adaptor_report.txt",
        logs = "results/{asmname}/5.quality_assessment/03.fcs/{asmname}/logs.jsonl",
        pipeline = "results/{asmname}/5.quality_assessment/03.fcs/{asmname}/pipeline_args.yaml",
        skipped = "results/{asmname}/5.quality_assessment/03.fcs/{asmname}/skipped_trims.jsonl",
        validate = "results/{asmname}/5.quality_assessment/03.fcs/{asmname}/validate_fasta.txt",
    log:
        "results/logs/5.quality_assessment/fcs-adaptor/{asmname}.log"
    benchmark:
        "results/benchmarks/5.quality_assessment/fcs-adaptor/{asmname}.txt"
    container:
        "docker://ncbi/fcs-adaptor:0.5.0"
    retries:
        3
    shell:
        "av_screen_x -o $(dirname {output.fcsadaptorrpt}) --euk {input} &> {log}"

rule visualise_fcs_adaptor:
    input:
        "results/{asmname}/5.quality_assessment/03.fcs/{asmname}/fcs_adaptor_report.txt"
    output:
        report("results/{asmname}/5.quality_assessment/03.fcs/{asmname}/fcs_adaptor_report.html",
            category="Contamination",
            caption="../../report/fcs-adaptor.rst",
            labels={"type": "fcs-adaptor", "assembly": "{asmname}"}),
    log:
        "results/logs/5.quality_assessment/visualise_fcs_adaptor/{asmname}.log"
    benchmark:
        "results/benchmarks/5.quality_assessment/visualise_fcs_adaptor/{asmname}.txt"
    conda:
        "../../envs/csvtotable.yaml"
    shell:
        "csvtotable -d $'\\t' {input} {output} &> {log}"
