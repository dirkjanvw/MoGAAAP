rule mash_distance:
    input:
        lambda wildcards: expand("final_output/{asmname}.full.fa", asmname=get_all_accessions_from_asmset(wildcards.asmset)),
    output:
        tsv = "results/{asmset}/3.quality_assessment/09.mash/{asmset}.tsv",
        csv = "results/{asmset}/3.quality_assessment/09.mash/{asmset}.csv",
    log:
        "results/logs/3.quality_assessment/mash_distance/{asmset}.log"
    benchmark:
        "results/benchmarks/3.quality_assessment/mash_distance/{asmset}.txt"
    threads:
        min(workflow.cores, 50)
    conda:
        "../../envs/ntsynt.yaml"
    shell:
        """
        (
        mash triangle -p {threads} {input} > {output.tsv}
        awk 'BEGIN{{FS = "\\t";OFS = ",";}} {{n=split($1,a,"/");$1=a[n];}} 1' {output.tsv} | awk 'BEGIN{{FS = OFS = ",";}} FNR!=1{{header=","$1""header; lines[$0];}} END{{print header; for (line in lines){{print line;}}}}' | sed 's/$/,/g' | sort > {output.csv}
        awk 'BEGIN{{FS = "\\t";OFS = ",";}} {{n=split($1,a,"/");$1=a[n];}} 1' {output.tsv} | awk 'BEGIN{{FS = OFS = ",";}} FNR!=1{{header[FNR]=$1; lines[FNR]=$0;}} END{{for (i in header){{printf ",%s",header[i];}} printf "\\n"; for (i in lines){{print lines[i];}}}}' | sed 's/$/,/g' > {output.csv}
        ) &> {log}
        """

rule mash_heatmap:
    input:
        "results/{asmset}/3.quality_assessment/09.mash/{asmset}.csv",
    output:
        report("results/{asmset}/3.quality_assessment/09.mash/{asmset}.pdf",
            category="Quality assessment",
            subcategory="Phylogeny",
            caption="../../report/mash.rst",
            labels={"type": "mash", "set": "{asmset}", "distance": "mash"}),
    log:
        "results/logs/3.quality_assessment/mash_heatmap/{asmset}.log"
    benchmark:
        "results/benchmarks/3.quality_assessment/mash_heatmap/{asmset}.txt"
    conda:
        "../../envs/rbase.yaml"
    shell:
        "Rscript workflow/scripts/pheatmap_lower_triangle.R {input} {output} &> {log}"
