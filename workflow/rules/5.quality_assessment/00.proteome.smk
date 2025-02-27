rule get_proteome:
    input:
        genome = "final_output/{asmname}.full.fa",
        annotation = "final_output/{asmname}.full.gff",
        config = "results/{asmname}/agat_config.yaml",
    output:
        "results/{asmname}/5.quality_assessment/proteome.pep.fa",
    log:
        "results/logs/5.quality_assessment/proteome/{asmname}.log",
    benchmark:
        "results/benchmarks/5.quality_assessment/proteome/{asmname}.txt",
    conda:
        "../../envs/agat.yaml"
    shell:
        "agat_sp_extract_sequences.pl -g {input.annotation} -f {input.genome} --cis --cfs -p -c {input.config} -o {output} &> {log}"

rule extract_isoforms:
    input:
        "results/{asmname}/5.quality_assessment/proteome.pep.fa",
    output:
        "results/{asmname}/5.quality_assessment/proteome.splice",
    log:
        "results/logs/5.quality_assessment/isoforms/{asmname}.log",
    benchmark:
        "results/benchmarks/5.quality_assessment/isoforms/{asmname}.txt",
    shell:
        "awk '/^>/{{split($2,a,\"=\"); i[a[2]] = i[a[2]]substr($1,2)\";\";}} END{{for (g in i){{print i[g];}}}}' {input} | sed 's/;$//g' > {output} 2> {log}"
