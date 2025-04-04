rule sans_prepare_genome:
    input:
        genomes = lambda wildcards: expand("final_output/{asmname}.full.fa", asmname=get_all_accessions_from_asmset(wildcards.asmset)),
    output:
        "results/{asmset}/3.quality_assessment/11.sans.list",
    log:
        "results/logs/3.quality_assessment/sans_prepare_genome/{asmset}.log"
    benchmark:
        "results/benchmarks/3.quality_assessment/sans_prepare_genome/{asmset}.txt"
    shell:
        """
        (
        for genome in {input.genomes}
        do
            name=$(basename ${{genome}} | rev | cut -d '.' -f 2- | rev)
            echo "${{name}} : ../../../${{genome}} ! 1" >> {output}
        done
        ) &> {log}
        """

rule sans:
    input:
        "results/{asmset}/3.quality_assessment/11.sans.list",
    output:
        splits = "results/{asmset}/3.quality_assessment/11.sans/{k}/{asmset}_b{bootstrap}.splits",
        bootstrap = "results/{asmset}/3.quality_assessment/11.sans/{k}/{asmset}_b{bootstrap}.splits.bootstrap",
        nexus = report("results/{asmset}/3.quality_assessment/11.sans/{k}/{asmset}_b{bootstrap}.nexus",
            category="Phylogeny",
            caption="../../report/sans.rst",
            labels={"type": "SANS", "set": "{asmset}", "k": "{k}", "bootstrap": "{bootstrap}"}),
    log:
        "results/logs/3.quality_assessment/sans/{k}/{asmset}_b{bootstrap}.log"
    benchmark:
        "results/benchmarks/3.quality_assessment/sans/{k}/{asmset}_b{bootstrap}.txt"
    threads:
        min(workflow.cores - 10, 10)
    container:
        "oras://ghcr.io/dirkjanvw/mogaaap/sans.v1.0.0:latest"
    shell:
        "SANS -i {input} -o {output.splits} -X {output.nexus} -f weakly -v -k {wildcards.k} -T {threads} -b {wildcards.bootstrap} &> {log}"
