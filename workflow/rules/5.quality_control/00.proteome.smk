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

rule extract_isoforms:
    input:
        "results/{asmname}/5.quality_control/proteome.pep.fa",
    output:
        "results/{asmname}/5.quality_control/proteome.splice",
    log:
        "results/logs/5.quality_control/isoforms/{asmname}.log",
    benchmark:
        "results/benchmark/5.quality_control/isoforms/{asmname}.txt",
    shell:
        "awk '/^>/{{split($2,a,\"=\"); i[a[2]] = i[a[2]]substr($1,2)\";\";}} END{{for (g in i){{print i[g];}}}}' {input} | sed 's/;$//g' > {output} 2> {log}"