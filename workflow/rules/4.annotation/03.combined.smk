rule combine_liftoff_helixer:
    input:
        liftoff = "results/{asmname}/4.annotation/01.liftoff/liftoff.gff_polished.fixed.valid_ORF.gff3", # fixed CDS phase and removed transcripts with invalid ORF
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

rule remove_attributes_gff:
    input:
        original = "results/{asmname}/4.annotation/03.combined/{asmname}.gff",
        config = "results/{asmname}/agat_config.yaml",
    output:
        temporary("results/{asmname}/4.annotation/03.combined/{asmname}.no_attr.gff"),
    log:
        "results/logs/4.annotation/remove_attributes_gff/{asmname}.log"
    benchmark:
        "results/benchmarks/4.annotation/remove_attributes_gff/{asmname}.txt"
    conda:
        "../../envs/agat.yaml"
    shell:
        "agat_sq_manage_attributes.pl --gff {input.original} --tag all_attributes --out {output} --config {input.config} &> {log}"

rule create_clean_gff:
    input:
        original = "results/{asmname}/4.annotation/03.combined/{asmname}.no_attr.gff",
        config = "results/{asmname}/agat_config.yaml",
    output:
        "results/{asmname}/4.annotation/03.combined/{asmname}.clean.gff",
    log:
        "results/logs/4.annotation/create_clean_gff/{asmname}.log"
    benchmark:
        "results/benchmarks/4.annotation/create_clean_gff/{asmname}.txt"
    params:
        species = lambda wildcards: get_species_name(wildcards),
    conda:
        "../../envs/agat.yaml"
    shell:
        "agat_sp_manage_IDs.pl --gff {input.original} --prefix {params.species} --tair --out {output} --config {input.config} &> {log}"
