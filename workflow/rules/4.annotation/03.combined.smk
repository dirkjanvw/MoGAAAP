rule gffcompare_liftoff_helixer:
    input:
        liftoff = "results/{asmname}/4.annotation/01.liftoff/liftoff.gff_polished",
        helixer = "results/{asmname}/4.annotation/02.helixer/helixer.gff",
    output:
        directory("results/{asmname}/4.annotation/03.combined/gffcompare_liftoff_helixer/"),
    log:
        "results/logs/4.annotation/gffcompare_liftoff_helixer/{asmname}.log"
    benchmark:
        "results/benchmarks/4.annotation/gffcompare_liftoff_helixer/{asmname}.txt"
    conda:
        "../../envs/gffcompare.yaml"
    shell:
        "gffcompare -r {input.liftoff} -o {output} {input.helixer} &> {log}"

rule create_helixer_kill_list:
    input:
        "results/{asmname}/4.annotation/03.combined/gffcompare_liftoff_helixer/gffcmp.tracking"
    output:
        "results/{asmname}/4.annotation/03.combined/helixer_kill_list.txt"
    log:
        "results/logs/4.annotation/create_helixer_kill_list/{asmname}.log"
    benchmark:
        "results/benchmarks/4.annotation/create_helixer_kill_list/{asmname}.txt"
    shell:  #only u, s, x, i, y, p types of overlap are kept
        "awk 'BEGIN{{FS = OFS = \"\t\";}} $4!~/^[usxiyp]/${{split($5,q,\"|\"); print q[2];}}' {input} > {output} 2> {log}"

rule filter_helixer:
    input:
        helixer = "results/{asmname}/4.annotation/02.helixer/helixer.gff",
        kill_list = "results/{asmname}/4.annotation/03.combined/helixer_kill_list.txt"
        config = "results/{asmname}/agat_config.yaml",
    output:
        "results/{asmname}/4.annotation/03.combined/helixer_filtered.gff"
    log:
        "results/logs/4.annotation/filter_helixer/{asmname}.log"
    benchmark:
        "results/benchmarks/4.annotation/filter_helixer/{asmname}.txt"
    conda:
        "../../envs/agat.yaml"
    shell:
        "agat_sp_filter_feature_from_kill_list.pl --gff {input.helixer} --kill_list {input.kill_list} --output {output} --config {input.config} &> {log}"

rule combine_liftoff_helixer:
    input:
        liftoff = "results/{asmname}/4.annotation/01.liftoff/liftoff.gff_polished",
        helixer = "results/{asmname}/4.annotation/03.combined/helixer_filtered.gff",
    output:
        "results/{asmname}/4.annotation/03.combined/{asmname}.gff"
    log:
        "results/logs/4.annotation/combine_liftoff_helixer/{asmname}.log"
    benchmark:
        "results/benchmarks/4.annotation/combine_liftoff_helixer/{asmname}.txt"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        "cat {input} | bedtools sort -i - > {output} 2> {log}"
