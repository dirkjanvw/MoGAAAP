include: "01.liftoff.smk"
include: "02.helixer.smk"
include: "03.combined.smk"

rule copy_annotation:
    input:
        full = lambda wildcards: "results/{asmname}/4.annotation/03.combined/{asmname}.gff"
            if PERFORM_ASSEMBLY
            else get_annotation_location(wildcards.asmname),
        coding = "results/{asmname}/4.annotation/03.combined/{asmname}.coding.gff",
        clean = "results/{asmname}/4.annotation/03.combined/{asmname}.clean.gff",
    output:
        full = "final_output/{asmname}.full.gff",
        coding = "final_output/{asmname}.full.coding.gff",
        clean = "final_output/{asmname}.full.clean.gff",
    log:
        "results/logs/4.annotation/copy_annotation/{asmname}.log"
    benchmark:
        "results/benchmarks/4.annotation/copy_annotation/{asmname}.txt"
    shell:
        """
        (
        cp {input.full} {output.full}
        cp {input.coding} {output.coding}
        cp {input.clean} {output.clean}
        ) &> {log}
        """

rule annotate:
    input:
        expand("final_output/{asmname}.full.gff", asmname=get_all_accessions()),
        expand("final_output/{asmname}.full.coding.gff", asmname=get_all_accessions()),
        expand("final_output/{asmname}.full.clean.gff", asmname=get_all_accessions()),
    output:
        touch("results/annotation.done")
