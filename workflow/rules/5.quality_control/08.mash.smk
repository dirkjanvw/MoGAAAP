rule mash_distance:
    input:
        lambda wildcards: expand("results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa", asmname=config["set"][wildcards.asmset]),
    output:
        tsv = "results/{asmset}/5.quality_control/mash/{asmset}.tsv",
        csv = "results/{asmset}/5.quality_control/mash/{asmset}.csv",
    log:
        "results/logs/5.quality_control/mash_distance/{asmset}.log"
    benchmark:
        "results/benchmarks/5.quality_control/mash_distance/{asmset}.txt"
    threads:
        min(workflow.cores, 50)
    conda:
        "../../envs/ntsynt.yaml"
    shell:
        """
        (
        mash triangle -p {threads} {input} > {output.tsv}
        sed 's/\\t/,/g' {output.tsv} | awk 'BEGIN{{FS = OFS = ",";}} FNR!=1{{header=","$1""header; lines[$0];}} END{{print header; for (line in lines){{print line;}}}}' | sed 's/$/,/g' | sort > {output.csv}
        sed 's/\\t/,/g' {output.tsv} | awk 'BEGIN{{FS = OFS = ",";}} FNR!=1{{header[FNR]=$1; lines[FNR]=$0;}} END{{for (i in header){{printf ",%s",header[i];}} printf "\\n"; for (i in lines){{print lines[i];}}}}' | sed 's/$/,/g' > {output.csv}
        ) &> {log}
        """

rule mash_heatmap:
    input:
        "results/{asmset}/5.quality_control/mash/{asmset}.csv",
    output:
        report("results/{asmset}/5.quality_control/08.mash/{asmset}.pdf", category="General", labels={"type": "mash", "set": "{asmset}", "distance": "mash"}),
    log:
        "results/logs/5.quality_control/mash_heatmap/{asmset}.log"
    benchmark:
        "results/benchmarks/5.quality_control/mash_heatmap/{asmset}.txt"
    conda:
        "../../envs/rbase.yaml"
    shell:
        "Rscript workflow/scripts/pheatmap_lower_triangle.R {input} {output} &> {log}"