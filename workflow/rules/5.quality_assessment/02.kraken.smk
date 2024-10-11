rule kraken2:
    input:
        "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa",
    output:
        out = "results/{asmname}/5.quality_assessment/02.kraken2/{asmname}.kraken2.out",
        report = "results/{asmname}/5.quality_assessment/02.kraken2/{asmname}.kraken2.report.txt",
    log:
        "results/logs/5.quality_assessment/kraken2/{asmname}.log"
    benchmark:
        "results/benchmarks/5.quality_assessment/kraken2/{asmname}.txt"
    params:
        db=config["kraken2_nt"]
    resources:
        gbmem=500  #should be just above what the kraken2 nt database actually uses
    threads:
        10
    conda:
        "../../envs/kraken2.yaml"
    shell:
        "kraken2 --db {params.db} --threads {threads} --output {output.out} --report {output.report} {input} > {output} 2> {log}"

rule krona:
    input:
        "results/{asmname}/5.quality_assessment/02.kraken2/{asmname}.kraken2.out",
    output:
        report("results/{asmname}/5.quality_assessment/02.kraken2/{asmname}.kraken2.krona.html",
            category="Contamination",
            caption="../../report/kraken.rst",
            labels={"type": "kraken2", "assembly": "{asmname}"}),
    log:
        "results/logs/5.quality_assessment/krona/{asmname}.log"
    benchmark:
        "results/benchmarks/5.quality_assessment/krona/{asmname}.txt"
    conda:
        "../../envs/kraken2.yaml"
    shell:
        """
        (
        output={output};
        krona=${{output%.html}};
        cut -f 2,3 {input} > $krona;
        ktImportTaxonomy $krona -o {output};
        ) &> {log}
        """