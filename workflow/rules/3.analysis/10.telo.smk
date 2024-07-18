rule identify_telomere:
    input:
        "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa",
    output:
        "results/{asmname}/3.analysis/10.telo/{asmname}.telo.bed",
    log:
        "results/logs/3.analysis/identify_telomere/{asmname}.log"
    benchmark:
        "results/benchmarks/3.analysis/identify_telomere/{asmname}.txt"
    params:
        telomere_motif = config["telomere_motif"]
    conda:
        "../../envs/seqtk.yaml"
    shell:
        "(seqtk telo -m {params.telomere_motif} {input} | sed '$ d' | cut -f -3 > {output}) 2> {log}"

rule visualise_telomere_locations:
    input:
        "results/{asmname}/3.analysis/10.telo/{asmname}.telo.bed",
    output:
        report("results/{asmname}/3.analysis/10.telo/{asmname}.telo.html",
            category="Analysis",
            caption="../../report/telo.rst",
            labels={"asmname": "{asmname}"}),
    log:
        "results/logs/3.analysis/visualise_telomere_locations/{asmname}.log"
    benchmark:
        "results/benchmarks/3.analysis/visualise_telomere_locations/{asmname}.txt"
    conda:
        "../../envs/csvtotable.yaml"
    shell:
        "csvtotable -d $'\\t' {input} {output} &> {log}"
