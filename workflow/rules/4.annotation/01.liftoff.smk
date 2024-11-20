rule liftoff:
    input:
        ref_annotation = lambda wildcards: config["reference_genomes"][get_reference_id(wildcards.asmname)]["annotation"],
        ref_genome = lambda wildcards: config["reference_genomes"][get_reference_id(wildcards.asmname)]["genome"],
        assembly = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa",
    output:
        ref_annotation = temporary("results/{asmname}/4.annotation/01.liftoff/reference.gff"),  #to prevent write permission issues
        ref_genome = temporary("results/{asmname}/4.annotation/01.liftoff/reference.fa"),  #to prevent write permission issues
        gff = "results/{asmname}/4.annotation/01.liftoff/liftoff.gff",
        polished = "results/{asmname}/4.annotation/01.liftoff/liftoff.gff_polished",
    log:
        "results/logs/4.annotation/liftoff/{asmname}.log"
    benchmark:
        "results/benchmarks/4.annotation/liftoff/{asmname}.txt"
    threads:
        10
    conda:
        "../../envs/liftoff.yaml"
    shell:
        """
        (
        ln -s $(realpath {input.ref_genome}) {output.ref_genome}
        ln -s $(realpath {input.ref_annotation}) {output.ref_annotation}
        liftoff -p {threads} -copies -cds -polish -u $(dirname {output.polished})/unmapped_features.txt -dir $(dirname {output.polished})/intermediate_files -o {output.gff} -g {output.ref_annotation} {input.assembly} {output.ref_genome}
        ) &> {log}
        """

rule fix_cds_phase_liftoff_gff_polished:
    input:
        polished = "results/{asmname}/4.annotation/01.liftoff/liftoff.gff_polished",
        assembly = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa",
        config = "results/{asmname}/agat_config.yaml",
    output:
        "results/{asmname}/4.annotation/01.liftoff/liftoff.gff_polished.fixed.gff",
    log:
        "results/logs/4.annotation/fixed_cds_phase_liftoff_gff_polished/{asmname}.log"
    benchmark:
        "results/benchmarks/4.annotation/fixed_cds_phase_liftoff_gff_polished/{asmname}.txt"
    conda:
        "../../envs/agat.yaml"
    shell:
        "agat_sp_fix_cds_phases.pl --gff {input.polished} --fasta {input.assembly} --output {output} --config {input.config} &> {log}"

rule filter_valid_ORF_mRNA:
    input:
        gff = "results/{asmname}/4.annotation/01.liftoff/liftoff.gff_polished.fixed.gff",
        config = "results/{asmname}/agat_config.yaml",
    output:
        gff = "results/{asmname}/4.annotation/01.liftoff/liftoff.gff_polished.fixed.valid_ORF.gff",
        report = "results/{asmname}/4.annotation/01.liftoff/liftoff.gff_polished.fixed.valid_ORF_report.txt",
        discarded = "results/{asmname}/4.annotation/01.liftoff/liftoff.gff_polished.fixed.valid_ORF_discarded.gff",
    log:
        "results/logs/4.annotation/filter_valid_ORF_mRNA/{asmname}.log"
    benchmark:
        "results/benchmarks/4.annotation/filter_valid_ORF_mRNA/{asmname}.txt"
    conda:
        "../../envs/agat.yaml"
    shell:
        "agat_sp_filter_feature_by_attribute_value.pl --gff {input.gff} -l level2 --attribute valid_ORF --value False --output {output.gff} --config {input.config} &> {log}"
