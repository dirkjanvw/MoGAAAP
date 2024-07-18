include: "01.liftoff.smk"
include: "02.helixer.smk"
include: "03.combined.smk"

rule copy_annotation:
    input:
        full = "results/{asmname}/4.annotation/03.combined/{asmname}.gff",
        coding = "results/{asmname}/4.annotation/03.combined/{asmname}.coding.gff",
    output:
        full = protected("final_output/{asmname}.full.gff"),
        coding = protected("final_output/{asmname}.coding.gff"),
    log:
        "results/logs/4.annotation/copy_annotation/{asmname}.log"
    benchmark:
        "results/benchmarks/4.annotation/copy_annotation/{asmname}.txt"
    shell:
        """
        (
        cp {input.full} {output.full}
        cp {input.coding} {output.coding}
        ) &> {log}
        """

rule annotate:
    input:
        expand("final_output/{asmname}.full.gff", asmname=get_all_accessions()),
    output:
        touch("results/annotation.done")
