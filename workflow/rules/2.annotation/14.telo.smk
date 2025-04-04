rule identify_telomere:
    input:
        "final_output/{asmname}.full.fa",
    output:
        "results/{asmname}/2.annotation/14.telo/{asmname}.telo.bed",
    log:
        "results/logs/2.annotation/identify_telomere/{asmname}.log"
    benchmark:
        "results/benchmarks/2.annotation/identify_telomere/{asmname}.txt"
    params:
        telomere_motif = config["telomere_motif"]
    conda:
        "../../envs/seqtk.yaml"
    shell:
        """
        (
        printf "sequence\\tstart\\tend\\n" > {output}
        seqtk telo -m {params.telomere_motif} {input} | sed '$ d' | cut -f -3 >> {output}
        ) 2> {log}
        """

rule visualise_telomere_locations:
    input:
        "results/{asmname}/2.annotation/14.telo/{asmname}.telo.bed",
    output:
        report("results/{asmname}/2.annotation/14.telo/{asmname}.telo.html",
            category="Custom annotation",
            caption="../../report/telo.rst",
            labels={"asmname": "{asmname}", "query_name": "telomere"}),
    log:
        "results/logs/2.annotation/visualise_telomere_locations/{asmname}.log"
    benchmark:
        "results/benchmarks/2.annotation/visualise_telomere_locations/{asmname}.txt"
    conda:
        "../../envs/csvtotable.yaml"
    shell:
        "csvtotable -d $'\\t' {input} {output} &> {log}"
