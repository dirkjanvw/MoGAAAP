def get_wgs_input(wildcards):
    if wildcards.sample in config["reads"]["hifi"][wildcards.asmname]:
        return config["reads"]["hifi"][wildcards.asmname][wildcards.sample]
    elif wildcards.sample in config["reads"]["ont"][wildcards.asmname]:
        return config["reads"]["ont"][wildcards.asmname][wildcards.sample]
    else:
        filelist = []
        for library in config["reads"]["illumina"][wildcards.asmname][wildcards.sample]:
            filelist.append(config["reads"]["illumina"][wildcards.asmname][wildcards.sample][library][0])
            filelist.append(config["reads"]["illumina"][wildcards.asmname][wildcards.sample][library][1])
        return filelist

rule meryl:
    input:
        get_wgs_input
    output:
        temporary(directory("results/{asmname}/5.quality_control/01.meryl_databases/{k}/{sample}.meryl")),  #relatively fast to compute and takes up a lot of space
    log:
        "results/logs/5.quality_control/meryl/{k}/{asmname}/{sample}.log"
    benchmark:
        "results/benchmarks/5.quality_control/meryl/{k}/{asmname}/{sample}.txt"
    threads:
        10
    conda:
        "../../envs/merqury.yaml"
    shell:
        "meryl count threads={threads} k={wildcards.k} output {output} {input} &> {log}"

rule merqury:
    input:
        meryl = "results/{asmname}/5.quality_control/01.meryl_databases/{k}/{sample}.meryl",
        genome = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa",
    output:
        temporary(directory("results/{asmname}/5.quality_control/01.merqury/{k}/{sample}/{asmname}.meryl")), #relatively fast to compute and takes up a lot of space
        bed = "results/{asmname}/5.quality_control/01.merqury/{k}/{sample}/{asmname}_only.bed",
        wig = "results/{asmname}/5.quality_control/01.merqury/{k}/{sample}/{asmname}_only.wig",
        stats = "results/{asmname}/5.quality_control/01.merqury/{k}/{sample}/{asmname}_vs_{sample}.completeness.stats",
        distonlyhist = "results/{asmname}/5.quality_control/01.merqury/{k}/{sample}/{asmname}_vs_{sample}.dist_only.hist",
        allqv = "results/{asmname}/5.quality_control/01.merqury/{k}/{sample}/{asmname}_vs_{sample}.qv",
        qv = "results/{asmname}/5.quality_control/01.merqury/{k}/{sample}/{asmname}_vs_{sample}.{sample}.qv",
        cnflplot = report("results/{asmname}/5.quality_control/01.merqury/{k}/{sample}/{asmname}_vs_{sample}.{sample}.spectra-cn.fl.png", category="K-mer completeness", labels={"type": "merqury", "assembly": "{asmname}", "wgs": "{sample}", "k": "{k}"}),
        cnhist = "results/{asmname}/5.quality_control/01.merqury/{k}/{sample}/{asmname}_vs_{sample}.{sample}.spectra-cn.hist",
        cnlnplot = "results/{asmname}/5.quality_control/01.merqury/{k}/{sample}/{asmname}_vs_{sample}.{sample}.spectra-cn.ln.png",
        cnstplot = "results/{asmname}/5.quality_control/01.merqury/{k}/{sample}/{asmname}_vs_{sample}.{sample}.spectra-cn.st.png",
        asmplot = "results/{asmname}/5.quality_control/01.merqury/{k}/{sample}/{asmname}_vs_{sample}.spectra-asm.fl.png",
        asmhist = "results/{asmname}/5.quality_control/01.merqury/{k}/{sample}/{asmname}_vs_{sample}.spectra-asm.hist",
        asmlnplot = "results/{asmname}/5.quality_control/01.merqury/{k}/{sample}/{asmname}_vs_{sample}.spectra-asm.ln.png",
        asmstplot = "results/{asmname}/5.quality_control/01.merqury/{k}/{sample}/{asmname}_vs_{sample}.spectra-asm.st.png",
        filt = "results/{asmname}/5.quality_control/01.merqury/{k}/{sample}/{asmname}.filt",
        hist = "results/{asmname}/5.quality_control/01.merqury/{k}/{sample}/{asmname}.hist",
        ploidy = "results/{asmname}/5.quality_control/01.merqury/{k}/{sample}/{asmname}.hist.ploidy",
    log:
        "results/logs/5.quality_control/01.merqury/{k}/{asmname}/{sample}.log"
    benchmark:
        "results/benchmarks/5.quality_control/01.merqury/{k}/{asmname}/{sample}.txt"
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
