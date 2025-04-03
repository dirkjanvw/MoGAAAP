def get_wgs_input(wildcards):
    if wildcards.wgstype == "hifi":
        return get_hifi(wildcards)
    elif wildcards.wgstype == "ont":
        return get_ont(wildcards)
    elif wildcards.wgstype == "illumina":
        return [get_illumina_1(wildcards), get_illumina_2(wildcards)]

rule meryl:
    input:
        get_wgs_input
    output:
        temporary(directory("results/{asmname}/5.quality_assessment/01.meryl_databases/{k}/{wgstype}.meryl")),  #relatively fast to compute and takes up a lot of space
    log:
        "results/logs/5.quality_assessment/meryl/{k}/{asmname}/{wgstype}.log"
    benchmark:
        "results/benchmarks/5.quality_assessment/meryl/{k}/{asmname}/{wgstype}.txt"
    threads:
        10
    conda:
        "../../envs/merqury.yaml"
    shell:
        "meryl count threads={threads} k={wildcards.k} output {output} {input} &> {log}"

rule merqury:
    input:
        meryl = "results/{asmname}/5.quality_assessment/01.meryl_databases/{k}/{wgstype}.meryl",
        genome = "final_output/{asmname}.full.fa",
    output:
        temporary(directory("results/{asmname}/5.quality_assessment/01.merqury/{k}/{wgstype}/{asmname}.full.meryl")), #relatively fast to compute and takes up a lot of space
        bed = "results/{asmname}/5.quality_assessment/01.merqury/{k}/{wgstype}/{asmname}.full_only.bed",
        wig = "results/{asmname}/5.quality_assessment/01.merqury/{k}/{wgstype}/{asmname}.full_only.wig",
        stats = "results/{asmname}/5.quality_assessment/01.merqury/{k}/{wgstype}/{asmname}_vs_{wgstype}.completeness.stats",
        distonlyhist = "results/{asmname}/5.quality_assessment/01.merqury/{k}/{wgstype}/{asmname}_vs_{wgstype}.dist_only.hist",
        allqv = "results/{asmname}/5.quality_assessment/01.merqury/{k}/{wgstype}/{asmname}_vs_{wgstype}.qv",
        qv = "results/{asmname}/5.quality_assessment/01.merqury/{k}/{wgstype}/{asmname}_vs_{wgstype}.{asmname}.full.qv",
        cnflplot = report("results/{asmname}/5.quality_assessment/01.merqury/{k}/{wgstype}/{asmname}_vs_{wgstype}.{asmname}.full.spectra-cn.fl.png",
            category="Quality assessment",
            subcategory="K-mer completeness",
            caption="../../report/merqury_plot.rst",
            labels={"type": "spectra-cn", "scope": "all sequences", "assembly": "{asmname}", "wgs": "{wgstype}", "k": "{k}"}),
        cnhist = "results/{asmname}/5.quality_assessment/01.merqury/{k}/{wgstype}/{asmname}_vs_{wgstype}.{asmname}.full.spectra-cn.hist",
        cnlnplot = "results/{asmname}/5.quality_assessment/01.merqury/{k}/{wgstype}/{asmname}_vs_{wgstype}.{asmname}.full.spectra-cn.ln.png",
        cnstplot = "results/{asmname}/5.quality_assessment/01.merqury/{k}/{wgstype}/{asmname}_vs_{wgstype}.{asmname}.full.spectra-cn.st.png",
        asmplot = "results/{asmname}/5.quality_assessment/01.merqury/{k}/{wgstype}/{asmname}_vs_{wgstype}.spectra-asm.fl.png",
        asmhist = "results/{asmname}/5.quality_assessment/01.merqury/{k}/{wgstype}/{asmname}_vs_{wgstype}.spectra-asm.hist",
        asmlnplot = "results/{asmname}/5.quality_assessment/01.merqury/{k}/{wgstype}/{asmname}_vs_{wgstype}.spectra-asm.ln.png",
        asmstplot = "results/{asmname}/5.quality_assessment/01.merqury/{k}/{wgstype}/{asmname}_vs_{wgstype}.spectra-asm.st.png",
        filt = "results/{asmname}/5.quality_assessment/01.merqury/{k}/{wgstype}/{wgstype}.filt",
        hist = "results/{asmname}/5.quality_assessment/01.merqury/{k}/{wgstype}/{wgstype}.hist",
        ploidy = "results/{asmname}/5.quality_assessment/01.merqury/{k}/{wgstype}/{wgstype}.hist.ploidy",
    log:
        "results/logs/5.quality_assessment/merqury/{k}/{asmname}/{wgstype}.log"
    benchmark:
        "results/benchmarks/5.quality_assessment/merqury/{k}/{asmname}/{wgstype}.txt"
    conda:
        "../../envs/merqury.yaml"
    shell:
        """
        (
        curr_dir=${{PWD}}
        prefix=$(basename {output.allqv} | rev | cut -d '.' -f 2- | rev)
        cd $(dirname {output.allqv})
        merqury.sh ${{curr_dir}}/{input.meryl} ${{curr_dir}}/{input.genome} ${{prefix}}
        ) &> {log}
        """

rule visualise_qv:
    input:
        "results/{asmname}/5.quality_assessment/01.merqury/{k}/{wgstype}/{asmname}_vs_{wgstype}.{asmname}.full.qv"
    output:
        tsv = "results/{asmname}/5.quality_assessment/01.merqury/{k}/{wgstype}/{asmname}_vs_{wgstype}.{asmname}.full.qv.tsv",
        html = report("results/{asmname}/5.quality_assessment/01.merqury/{k}/{wgstype}/{asmname}_vs_{wgstype}.{asmname}.full.qv.html",
            category="Quality assessment",
            subcategory="K-mer completeness",
            caption="../../report/merqury_qv.rst",
            labels={"type": "QV", "scope": "per sequence", "assembly": "{asmname}",
                    "wgs": "{wgstype}", "k": "{k}"}),
    log:
        "results/logs/5.quality_assessment/visualise_qv/{k}/{asmname}/{wgstype}.log"
    benchmark:
        "results/benchmarks/5.quality_assessment/visualise_qv/{k}/{asmname}/{wgstype}.txt"
    conda:
        "../../envs/csvtotable.yaml"
    shell:
        """
        (
        awk 'BEGIN{{FS = OFS = "\\t"; print "Sequence", "K-mers unique for assembly", "K-mers in both assembly and read set", "QV", "Error rate";}} {{print;}}' {input} > {output.tsv}
        csvtotable -d $'\\t' {output.tsv} {output.html}
        ) &> {log}
        """
