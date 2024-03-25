rule fcs_gx:
    input:
        "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa"
    output:
    report = report("results/{asmname}/5.quality_control/03.fcs/{asmname}.fcs_gx_report.txt", category="Contamination", labels={"type": "fcs-gx", "assembly": "{asmname}"}),
    taxonomy = "results/{asmname}/5.quality_control/03.fcs/{asmname}.taxonomy.rpt",
    log:
        "results/logs/5.quality_control/fcs-gx/{asmname}.log"
    benchmark:
        "results/benchmarks/5.quality_control/fcs-gx/{asmname}.txt"
    params:
        taxid = lambda wildcards: config["taxid"][wildcards.asmname],
        gxdb = config["gxdb"],
    resources:
        gbmem = 500
    container:
        "docker://ncbi/fcs-gx:0.5.0"
    shell:
        "run_gx --fasta {input} --tax-id {params.taxid} --gx-db {params.gxdb} --out-basename {wildcards.asmname} --out-dir $(dirname {output.taxonomy}) &> {log}"

rule fcs_adaptor:
    input:
        "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa"
    output:
        cleaned = "results/{asmname}/5.quality_control/03.fcs/{asmname}/cleaned_sequences/{asmname}.fa",
        calls = "results/{asmname}/5.quality_control/03.fcs/{asmname}/combined.calls.jsonl",
        fcslog = "results/{asmname}/5.quality_control/03.fcs/{asmname}/fcs.log",
        fcsadaptorlog = "results/{asmname}/5.quality_control/03.fcs/{asmname}/fcs_adaptor.log",
        fcsadaptorrpt = report("results/{asmname}/5.quality_control/03.fcs/{asmname}/fcs_adaptor_report.txt", category="Contamination", labels={"type": "fcs-adaptor", "assembly": "{asmname}"}),
        logs = "results/{asmname}/5.quality_control/03.fcs/{asmname}/logs.jsonl",
        pipeline = "results/{asmname}/5.quality_control/03.fcs/{asmname}/pipeline_args.yaml",
        skipped = "results/{asmname}/5.quality_control/03.fcs/{asmname}/skipped_trims.jsonl",
        validate = "results/{asmname}/5.quality_control/03.fcs/{asmname}/validate_fasta.txt",
    log:
        "results/logs/5.quality_control/fcs-adaptor/{asmname}.log"
    benchmark:
        "results/benchmarks/5.quality_control/fcs-adaptor/{asmname}.txt"
    container:
        "docker://ncbi/fcs-adaptor:0.5.0"
    shell:
        "av_screen_x -o $(dirname {output.fcsadaptorrpt}) --euk {input} &> {log}"
