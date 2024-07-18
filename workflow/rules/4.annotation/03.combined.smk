rule combine_liftoff_helixer:
    input:
        liftoff = "results/{asmname}/4.annotation/01.liftoff/liftoff.gff_polished",
        helixer = "results/{asmname}/4.annotation/02.helixer/helixer.gff",
        config = "results/{asmname}/agat_config.yaml",
    output:
        "results/{asmname}/4.annotation/03.combined/{asmname}.gff",
    log:
        "results/logs/4.annotation/combine_liftoff_helixer/{asmname}.log"
    benchmark:
        "results/benchmarks/4.annotation/combine_liftoff_helixer/{asmname}.txt"
    conda:
        "../../envs/agat.yaml"
    shell:
        "agat_sp_complement_annotations.pl --ref {input.liftoff} --add {input.helixer} --out {output} --config {input.config} &> {log}"

rule identify_coding_genes:
    input:
        gff = "results/{asmname}/4.annotation/03.combined/{asmname}.gff",
        config = "results/{asmname}/agat_config.yaml",
    output:
        temp = temporary("results/{asmname}/4.annotation/03.combined/{asmname}.coding.list_ID"),
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
