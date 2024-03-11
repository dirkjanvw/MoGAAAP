rule combine_liftoff_helixer:
    input:
        liftoff = "results/{asmname}/4.annotation/01.liftoff/liftoff.gff_polished",
        helixer = "results/{asmname}/4.annotation/02.helixer/helixer.gff",
        config = "results/{asmname}/agat_config.yaml",
    output:
        protected("results/{asmname}/4.annotation/03.combined/{asmname}.gff"),
    log:
        "results/logs/4.annotation/combine_liftoff_helixer/{asmname}.log"
    benchmark:
        "results/benchmarks/4.annotation/combine_liftoff_helixer/{asmname}.txt"
    conda:
        "../../envs/agat.yaml"
    shell:
        "agat_sp_complement_annotations.pl --ref {input.liftoff} --add {input.helixer} --out {output} --config {input.config} &> {log}"
