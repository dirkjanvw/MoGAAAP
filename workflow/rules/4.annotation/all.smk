include: "01.liftoff.smk"
include: "02.helixer.smk"
include: "03.combined.smk"

rule link_annotation:
    input:
        "results/{asmname}/4.annotation/03.combined/{asmname}.gff"
    output:
        "final_output/{asmname}.full.gff"
    log:
        "results/logs/4.annotation/link_annotation/{asmname}.log"
    benchmark:
        "results/benchmarks/4.annotation/link_annotation/{asmname}.txt"
    shell:
        "cp $(realpath {input}) {output} &> {log}"

rule annotate:
    input:
        expand("final_output/{asmname}.full.gff", asmname=get_all_accessions()),
    output:
        touch("results/annotation.done")
