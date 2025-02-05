rule sans_prepare_genome:
    input:
        genomes = lambda wildcards: expand("results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa", asmname=get_all_accessions_from_asmset(wildcards.asmset)),
    output:
        "results/{asmset}/5.quality_assessment/11.sans.list",
    log:
        "results/logs/5.quality_assessment/sans_prepare_genome/{asmset}.log"
    benchmark:
        "results/benchmarks/5.quality_assessment/sans_prepare_genome/{asmset}.txt"
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
        "results/{asmset}/5.quality_assessment/11.sans.list",
    output:
        splits = "results/{asmset}/5.quality_assessment/11.sans/{k}/{asmset}_b{bootstrap}.splits",
        bootstrap = "results/{asmset}/5.quality_assessment/11.sans/{k}/{asmset}_b{bootstrap}.splits.bootstrap",
    log:
        "results/logs/5.quality_assessment/sans/{k}/{asmset}_b{bootstrap}.log"
    benchmark:
        "results/benchmarks/5.quality_assessment/sans/{k}/{asmset}_b{bootstrap}.txt"
    threads:
        min(workflow.cores - 10, 10)
    container:
        "oras://ghcr.io/dirkjanvw/mogaaap/sans.968e2d35be2a9f7fca75f664b004ef9cb32dd3e0:latest"
    shell:
        "SANS -i {input} -o {output.splits} -f weakly -v -k {wildcards.k} -T {threads} -b {wildcards.bootstrap} &> {log}"

rule sans_to_nexus:
    input:
        filelist = "results/{asmset}/5.quality_assessment/11.sans.list",
        splits = "results/{asmset}/5.quality_assessment/11.sans/{k}/{asmset}_b{bootstrap}.splits",
        bootstrap = "results/{asmset}/5.quality_assessment/11.sans/{k}/{asmset}_b{bootstrap}.splits.bootstrap",
    output:
        report("results/{asmset}/5.quality_assessment/11.sans/{k}/{asmset}_b{bootstrap}.nexus",
            category="Phylogeny",
            caption="../../report/sans.rst",
            labels={"type": "SANS", "set": "{asmset}", "k": "{k}", "bootstrap": "{bootstrap}"}),
    log:
        "results/logs/5.quality_assessment/sans/{k}/{asmset}_b{bootstrap}.nexus.log"
    benchmark:
        "results/benchmarks/5.quality_assessment/sans/{k}/{asmset}_b{bootstrap}.nexus.txt"
    container:
        "oras://ghcr.io/dirkjanvw/mogaaap/sans.968e2d35be2a9f7fca75f664b004ef9cb32dd3e0:latest"
    shell:
        "sans2conf_nexus.py {input.splits} {input.bootstrap} {input.filelist} > {output} 2> {log}"
