rule filter_sequences:
    input:
        "results/{asmname}/assembly/01.hifiasm/{asmname}.fa"
    output:
        "results/{asmname}/assembly/01.hifiasm/{asmname}.min{minlen}.fa"
    log:
        "results/logs/04.contigs_filtered/filter_sequences/{asmname}.{minlen}.log"
    benchmark:
        "results/benchmarks/04.contigs_filtered/filter_sequences/{asmname}.{minlen}.txt"
    conda:
        "../../envs/bioawk.yaml"
    shell:
        "bioawk -c fastx '{{ if(length($seq) >= {wildcards.minlen}) {{print \">\"$name; print $seq }} }}' {input} > {output} 2> {log}"

rule sort_sequences:
    input:
        "results/{asmname}/assembly/01.hifiasm/{asmname}.min{minlen}.fa"
    output:
        "results/{asmname}/assembly/01.hifiasm/{asmname}.min{minlen}.sorted.fa"
    log:
        "results/logs/04.contigs_filtered/sort_sequences/{asmname}.log"
    benchmark:
        "results/benchmarks/04.contigs_filtered/sort_sequences/{asmname}.txt"
    conda:
        "../../envs/bioawk.yaml"
    shell:
        """
        (
        bioawk -c fastx '{{ print ">"$name, length($seq), $seq }}' {input} | \
        sort -k2nr | \
        perl -pe 's/\\t/  L:/' | \
        perl -pe 's/\\t/ \\n/' > {output}
        ) &> {log}
        """

rule add_prefix:
    input:
        "results/{asmname}/assembly/01.hifiasm/{asmname}.min{minlen}.sorted.fa"
    output:
        "results/{asmname}/assembly/01.hifiasm/{asmname}.min{minlen}.sorted.renamed.fa"
    log:
        "results/logs/05.contigs_renamed/add_prefix/{asmname}.{minlen}.log"
    benchmark:
        "results/benchmarks/05.contigs_renamed/add_prefix/{asmname}.{minlen}.txt"
    params:
        prefix=lambda wildcards: f"{wildcards.asmname}_C"
    conda:
        "../../envs/bioawk.yaml"
    shell:
        """
        (
        bioawk -c fastx '{{ print ">{params.prefix}" $name; print $seq }}' {input} | \
        perl -pe 's/ptg00//' > {output}
        ) &> {log}
        """

