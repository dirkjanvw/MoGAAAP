rule pangrowth_hist:
    input:
        lambda wildcards: expand("final_output/{asmname}.full.fa", asmname=get_all_accessions_from_asmset(wildcards.asmset)),
    output:
        "results/{asmset}/3.quality_assessment/12.pangrowth/{k}/hist.txt",
    log:
        "results/logs/3.quality_assessment/pangrowth_hist/{k}/{asmset}.log"
    benchmark:
        "results/benchmarks/3.quality_assessment/pangrowth_hist/{k}/{asmset}.txt"
    threads:
        min(workflow.cores - 10, 10)
    container:
        "oras://ghcr.io/dirkjanvw/mogaaap/pangrowth.71d67bde89326644f6718c82ec2ee7b751f3080b:latest"
    shell:
        "pangrowth hist -k {wildcards.k} -t {threads} {input} > {output} 2> {log}"

rule pangrowth_hist_plot:
    input:
        "results/{asmset}/3.quality_assessment/12.pangrowth/{k}/hist.txt",
    output:
        report("results/{asmset}/3.quality_assessment/12.pangrowth/{k}/hist.pdf",
            category="Quality assessment",
            subcategory="Pangrowth",
            caption="../../report/pangrowth_hist.rst",
            labels={"type": "histogram", "set": "{asmset}", "k": "{k}"}),
    log:
        "results/logs/3.quality_assessment/pangrowth_hist_plot/{asmset}/{k}.log"
    benchmark:
        "results/benchmarks/3.quality_assessment/pangrowth_hist_plot/{asmset}/{k}.txt"
    container:
        "oras://ghcr.io/dirkjanvw/mogaaap/pangrowth.71d67bde89326644f6718c82ec2ee7b751f3080b:latest"
    shell:
        "plot_single_hist.py {input} {output} &> {log}"

rule pangrowth_growth:
    input:
        "results/{asmset}/3.quality_assessment/12.pangrowth/{k}/hist.txt",
    output:
        "results/{asmset}/3.quality_assessment/12.pangrowth/{k}/growth.txt",
    log:
        "results/logs/3.quality_assessment/pangrowth_growth/{asmset}/{k}.log"
    benchmark:
        "results/benchmarks/3.quality_assessment/pangrowth_growth/{asmset}/{k}.txt"
    container:
        "oras://ghcr.io/dirkjanvw/mogaaap/pangrowth.71d67bde89326644f6718c82ec2ee7b751f3080b:latest"
    shell:
        "pangrowth growth -h {input} > {output} 2> {log}"

rule pangrowth_growth_plot:
    input:
        "results/{asmset}/3.quality_assessment/12.pangrowth/{k}/growth.txt",
    output:
        report("results/{asmset}/3.quality_assessment/12.pangrowth/{k}/growth.pdf",
            category="Quality assessment",
            subcategory="Pangrowth",
            caption="../../report/pangrowth_growth.rst",
            labels={"type": "growth", "set": "{asmset}", "k": "{k}"}),
    log:
        "results/logs/3.quality_assessment/pangrowth_growth_plot/{asmset}/{k}.log"
    benchmark:
        "results/benchmarks/3.quality_assessment/pangrowth_growth_plot/{asmset}/{k}.txt"
    container:
        "oras://ghcr.io/dirkjanvw/mogaaap/pangrowth.71d67bde89326644f6718c82ec2ee7b751f3080b:latest"
    shell:
        "plot_growth.py {input} {output} &> {log}"

rule pangrowth_core:
    input:
        "results/{asmset}/3.quality_assessment/12.pangrowth/{k}/hist.txt",
    output:
        "results/{asmset}/3.quality_assessment/12.pangrowth/{k}/core.txt",
    log:
        "results/logs/3.quality_assessment/pangrowth_core/{asmset}/{k}.log"
    benchmark:
        "results/benchmarks/3.quality_assessment/pangrowth_core/{asmset}/{k}.txt"
    container:
        "oras://ghcr.io/dirkjanvw/mogaaap/pangrowth.71d67bde89326644f6718c82ec2ee7b751f3080b:latest"
    shell:
        "pangrowth core -h {input} > {output} 2> {log}"

rule pangrowth_core_plot:
    input:
        "results/{asmset}/3.quality_assessment/12.pangrowth/{k}/core.txt",
    output:
        report("results/{asmset}/3.quality_assessment/12.pangrowth/{k}/core.pdf",
            category="Quality assessment",
            subcategory="Pangrowth",
            caption="../../report/pangrowth_core.rst",
            labels={"type": "core", "set": "{asmset}", "k": "{k}"}),
    log:
        "results/logs/3.quality_assessment/pangrowth_core_plot/{asmset}/{k}.log"
    benchmark:
        "results/benchmarks/3.quality_assessment/pangrowth_core_plot/{asmset}/{k}.txt"
    container:
        "oras://ghcr.io/dirkjanvw/mogaaap/pangrowth.71d67bde89326644f6718c82ec2ee7b751f3080b:latest"
    shell:
        "plot_core.py {input} {output} &> {log}"
