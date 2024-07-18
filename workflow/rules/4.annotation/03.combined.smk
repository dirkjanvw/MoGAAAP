rule gffcompare_liftoff_helixer:
    input:
        liftoff = "results/{asmname}/4.annotation/01.liftoff/liftoff.gff_polished",
        helixer = "results/{asmname}/4.annotation/02.helixer/helixer.gff",
    output:
        "results/{asmname}/4.annotation/03.combined/gffcompare_liftoff_helixer/gffcmp.tracking",
    log:
        "results/logs/4.annotation/gffcompare_liftoff_helixer/{asmname}.log"
    benchmark:
        "results/benchmarks/4.annotation/gffcompare_liftoff_helixer/{asmname}.txt"
    conda:
        "../../envs/gffcompare.yaml"
    shell:
        "gffcompare -r {input.liftoff} -o $(dirname {output}) {input.helixer} &> {log}"

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
        kill_list = "results/{asmname}/4.annotation/03.combined/helixer_kill_list.txt",
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

rule identify_coding_genes:
    input:
        gff = "results/{asmname}/4.annotation/03.combined/{asmname}.gff",
        config = "results/{asmname}/agat_config.yaml",
    output:
        temp = temporary("results/{asmname}/4.annotation/03.combined/{asmname}.coding_ID.list"),
        list = "results/{asmname}/4.annotation/03.combined/{asmname}.coding.list",
    log:
        "results/logs/4.annotation/identify_coding_genes/{asmname}.log"
    benchmark:
        "results/benchmarks/4.annotation/identify_coding_genes/{asmname}.txt"
    conda:
        "../../envs/agat.yaml"
    shell:
        """
        (
        agat_sp_extract_attributes.pl --gff {input.gff} -p mRNA --attribute ID --output {output.list} --config {input.config}
        cut -d ',' -f 1 {output.temp} > {output.list}
        ) &> {log}
        """

rule filter_coding_genes:
    input:
        original = "results/{asmname}/4.annotation/03.combined/{asmname}.gff",
        coding_genes = "results/{asmname}/4.annotation/03.combined/{asmname}.coding.list",
        config = "results/{asmname}/agat_config.yaml",
    output:
        "results/{asmname}/4.annotation/03.combined/{asmname}.coding.gff",
    log:
        "results/logs/4.annotation/filter_coding_genes/{asmname}.log"
    benchmark:
        "results/benchmarks/4.annotation/filter_coding_genes/{asmname}.txt"
    conda:
        "../../envs/agat.yaml"
    shell:
        "agat_sp_filter_feature_from_keep_list.pl --gff {input.original} --keep_list {input.coding_genes} --output {output} --config {input.config} &> {log}"
