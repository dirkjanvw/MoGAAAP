rule list_pantools_proteomes:
    input:
        expand("results/{asmname}/5.quality_control/proteome.pep.fa", asmname=get_all_accessions()),
    output:
        "results/{asmset}/5.quality_control/04.pantools/panproteome.list"
    log:
        "results/logs/5.quality_control/pantools/list_proteomes/{asmset}.log"
    benchmark:
        "results/benchmarks/5.quality_control/pantools/list_proteomes/{asmset}.txt"
    shell:
        "ls {input} > {output} 2> {log}"

rule table_pantools_proteomes:
    input:
        "results/{asmset}/5.quality_control/04.pantools/panproteome.list"
    output:
        "results/{asmset}/5.quality_control/04.pantools/panproteome.tsv"
    log:
        "results/logs/5.quality_control/pantools/table_proteomes/{asmset}.log"
    benchmark:
        "results/benchmarks/5.quality_control/pantools/table_proteomes/{asmset}.txt"
    shell:
        "awk 'BEGIN{{OFS = \"\\t\";}} {{n=split($1,a,\"/\"); species=a[n]; sub(/\\.pep\\.fa/, \"\", species); print species,FNR;}}' {input} > {output} 2> {log}"

rule panproteome_build:
    input:
        "results/{asmset}/5.quality_control/04.pantools/panproteome.list"
    output:
        directory("results/{asmset}/5.quality_control/04.pantools/panproteome_DB/databases")
    log:
        "results/logs/5.quality_control/pantools/build_panproteome/{asmset}.log"
    benchmark:
        "results/benchmarks/5.quality_control/pantools/build_panproteome/{asmset}.txt"
    params:
        tmpdir = config["tmpdir"],
        jvm = config["jvm"],
    resources:
        pantools = 1  #to prevent multiple pantools runs at the same time
    conda:
        "../../envs/pantools.yaml"
    shell:
        """
        (
        mkdir -p {params.tmpdir}
        pantools {params.jvm} build_panproteome --force {params.tmpdir}/panproteome_DB {input}
        mv {params.tmpdir}/panproteome_DB $(dirname $(dirname {output}))
        ) &> {log}"""

rule panproteome_group:
    input:
        rules.panproteome_build.output
    output:
        "results/{asmset}/5.quality_control/04.pantools/panproteome_groups_DB/pantools_homology_groups.txt"
    log:
        "results/logs/5.quality_control/pantools/group/{asmset}.log"
    benchmark:
        "results/benchmarks/5.quality_control/pantools/group/{asmset}.txt"
    params:
        grouping = config["pantools_grouping"],
        tmpdir = config["tmpdir"],
        jvm = config["jvm"],
    threads:
        max(50, workflow.cores // 2 - 1)
    resources:
        pantools = 1  #to prevent multiple pantools runs at the same time
    conda:
        "../../envs/pantools.yaml"
    shell:
        """
        (
        mkdir -p {params.tmpdir}
        cp -vr $(dirname {input}) {params.tmpdir}/panproteome_groups_DB
        pantools {params.jvm} group --relaxation {params.grouping} --threads {threads} {params.tmpdir}/panproteome_groups_DB
        mv {params.tmpdir}/panproteome_groups_DB $(dirname $(dirname {output}))/
        ) &> {log}"""

rule panproteome_gene_classification:
    input:
        "results/{asmset}/5.quality_control/04.pantools/panproteome_groups_DB/pantools_homology_groups.txt"
    output:
        groups = "results/{asmset}/5.quality_control/04.pantools/panproteome_groups_DB/gene_classification/classified_groups.csv",
        rscript = "results/{asmset}/5.quality_control/04.pantools/panproteome_groups_DB/gene_classification/upset/upset_plot.R", #will only be an Rscript if less than 11 genomes
    log:
        "results/logs/5.quality_control/pantools/gene_classification/{asmset}.log"
    benchmark:
        "results/benchmarks/5.quality_control/pantools/gene_classification/{asmset}.txt"
    params:
        jvm = config["jvm"],
    resources:
        pantools = 1  #to prevent multiple pantools runs at the same time
    conda:
        "../../envs/pantools.yaml"
    shell:
        """
        (
        pantools {params.jvm} gene_classification $(dirname {input})
        [ ! -f {output.rscript} ] && touch {output.rscript}
        ) &> {log}
        """

rule panproteome_plot_upset:
    input:
        "results/{asmset}/5.quality_control/04.pantools/panproteome_groups_DB/gene_classification/upset/upset_plot.R",
    output:
        report("results/{asmset}/5.quality_control/04.pantools/panproteome_groups_DB/gene_classification/upset/output/genomes.pdf", category="PanTools", labels={"type": "upset plot", "set": "{asmset}"}),
    log:
        "results/logs/5.quality_control/pantools/upset_plot/{asmset}.log"
    benchmark:
        "results/benchmarks/5.quality_control/pantools/upset_plot/{asmset}.txt"
    conda:
        "../../envs/rbase.yaml"  #cannot use pantools environment due to conflicts?
    shell:
        "Rscript {input} &> {log}"

rule panproteome_pangenome_structure:
    input:
        hm = "results/{asmset}/5.quality_control/04.pantools/panproteome_groups_DB/pantools_homology_groups.txt",
        gc = "results/{asmset}/5.quality_control/04.pantools/panproteome_groups_DB/gene_classification/classified_groups.csv",  #put this here to prevent parallel pantools execution
    output:
        "results/{asmset}/5.quality_control/04.pantools/panproteome_groups_DB/pangenome_size/gene/pangenome_growth.R",
    log:
        "results/logs/5.quality_control/pantools/pangenome_structure/{asmset}.log"
    benchmark:
        "results/benchmarks/5.quality_control/pantools/pangenome_structure/{asmset}.txt"
    params:
        jvm = config["jvm"],
    threads:
        50
    resources:
        pantools = 1  #to prevent multiple pantools runs at the same time
    conda:
        "../../envs/pantools.yaml"
    shell:
        "pantools {params.jvm} pangenome_structure --threads {threads} $(dirname {input.hm}) &> {log}"

rule panproteome_plot_pangenome_growth:
    input:
        "results/{asmset}/5.quality_control/04.pantools/panproteome_groups_DB/pangenome_size/gene/pangenome_growth.R",
    output:
        report("results/{asmset}/5.quality_control/04.pantools/panproteome_groups_DB/pangenome_size/gene/core_accessory_unique_growth.png", category="PanTools", labels={"type": "growth (core, accessory, unique)", "set": "{asmset}"}),
        report("results/{asmset}/5.quality_control/04.pantools/panproteome_groups_DB/pangenome_size/gene/core_dispensable_growth.png", category="PanTools", labels={"type": "growth (core, dispensable)", "set": "{asmset}"}),
        report("results/{asmset}/5.quality_control/04.pantools/panproteome_groups_DB/pangenome_size/gene/core_dispensable_total_growth.png", category="PanTools", labels={"type": "growth (core, dispensable, total)", "set": "{asmset}"}),
    log:
        "results/logs/5.quality_control/pantools/pangenome_growth/{asmset}.log"
    benchmark:
        "results/benchmarks/5.quality_control/pantools/pangenome_growth/{asmset}.txt"
    params:
        jvm = config["jvm"],
    conda:
        "../../envs/pantools.yaml"
    shell:
        "Rscript {input} &> {log}"