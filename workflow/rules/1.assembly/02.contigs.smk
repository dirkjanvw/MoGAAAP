rule filter_sequences:
    input:
        "results/{asmname}/1.assembly/01.hifiasm/{asmname}.fa"
    output:
        "results/{asmname}/1.assembly/01.hifiasm/{asmname}.min{minlen}.fa"
    log:
        "results/logs/1.assembly/filter_sequences/{asmname}.{minlen}.log"
    benchmark:
        "results/benchmarks/1.assembly/filter_sequences/{asmname}.{minlen}.txt"
    conda:
        "../../envs/bioawk.yaml"
    shell:
        "bioawk -c fastx '{{ if(length($seq) >= {wildcards.minlen}) {{print \">\"$name; print $seq }} }}' {input} > {output} 2> {log}"

rule sort_sequences:
    input:
        "results/{asmname}/1.assembly/01.hifiasm/{asmname}.min{minlen}.fa"
    output:
        "results/{asmname}/1.assembly/01.hifiasm/{asmname}.min{minlen}.sorted.fa"
    log:
        "results/logs/1.assembly/sort_sequences/{asmname}.min{minlen}.log"
    benchmark:
        "results/benchmarks/1.assembly/sort_sequences/{asmname}.min{minlen}.txt"
    conda:
        "../../envs/bioawk.yaml"
    shell:
        "bioawk -c fastx '{{ print \">\"$name, length($seq), $seq }}' {input} | sort -k2nr | perl -pe 's/\\t/  L:/' | perl -pe 's/\\t/ \\n/' > {output} 2> {log}"

rule add_prefix:
    input:
        "results/{asmname}/1.assembly/01.hifiasm/{asmname}.min{minlen}.sorted.fa"
    output:
        "results/{asmname}/1.assembly/01.hifiasm/{asmname}.min{minlen}.sorted.renamed.fa"
    log:
        "results/logs/1.assembly/add_prefix/{asmname}.min{minlen}.log"
    benchmark:
        "results/benchmarks/1.assembly/add_prefix/{asmname}.min{minlen}.txt"
    params:
        prefix=lambda wildcards: f"{wildcards.asmname}_C"
    conda:
        "../../envs/bioawk.yaml"
    shell:
        "bioawk -c fastx '{{ print \">{params.prefix}\" $name; print $seq }}' {input} | perl -pe 's/ptg00//' > {output} 2> {log}"
