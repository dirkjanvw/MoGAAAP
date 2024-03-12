rule get_proteome:
    input:
        genome = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa",
        annotation = "results/{asmname}/4.annotation/03.combined/{asmname}.gff",
        config = "results/{asmname}/agat_config.yaml",
    output:
        "results/{asmname}/5.quality_control/proteome.pep.fa",
    log:
        "results/logs/5.quality_control/proteome/{asmname}.log",
    benchmark:
        "results/benchmark/5.quality_control/proteome/{asmname}.txt",
    conda:
        "../../envs/agat.yaml"
    shell:
        "agat_sp_extract_sequences.pl -g {input.annotation} -f {input.genome} --cis --cfs -p -c {input.config} -o {output} &> {log}"