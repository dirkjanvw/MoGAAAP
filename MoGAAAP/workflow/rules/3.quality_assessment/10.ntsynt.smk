rule ntsynt:
    input:
        genomes = lambda wildcards: expand("final_output/{asmname}.full.fa", asmname=get_all_accessions_from_asmset(wildcards.asmset)),
        divergence = "results/{asmset}/3.quality_assessment/09.mash/{asmset}.tsv",
    output:
        blocks = "results/{asmset}/3.quality_assessment/10.ntsynt/{asmset}.k{mink}.w{minw}.synteny_blocks.tsv",
        commonbf = "results/{asmset}/3.quality_assessment/10.ntsynt/{asmset}.k{mink}.w{minw}.common.bf",
        mxdot = "results/{asmset}/3.quality_assessment/10.ntsynt/{asmset}.k{mink}.w{minw}.mx.dot",
        precolblocks = "results/{asmset}/3.quality_assessment/10.ntsynt/{asmset}.k{mink}.w{minw}.pre-collinear-merge.synteny_blocks.tsv",
    log:
        "results/logs/3.quality_assessment/ntsynt/{asmset}.k{mink}.w{minw}.log"
    benchmark:
        "results/benchmarks/3.quality_assessment/ntsynt/{asmset}.k{mink}.w{minw}.txt"
    threads:
        min(workflow.cores, 50)
    conda:
        "../../envs/ntsynt.yaml"
    shell:
        """
        (
        divergence=$(awk '{{for (i=2;i<=NF;i++){{if ($i>d){{d=$i;}}}}}} END{{printf "%.2f",d*100;}}' {input.divergence})
        ln -sf $(realpath {input.genomes}) $(dirname {output.blocks})/
        cd $(dirname {output.blocks})
        ntSynt -p {wildcards.asmset}.k{wildcards.mink}.w{wildcards.minw} -k {wildcards.mink} -w {wildcards.minw} -d ${{divergence}} -f -t {threads} $(for file in {input.genomes}; do basename ${{file}}; done)
        ) &> {log}
        """

rule visualise_ntsynt:
    input:
        blocks = "results/{asmset}/3.quality_assessment/10.ntsynt/{asmset}.k{mink}.w{minw}.synteny_blocks.tsv",
        fais = lambda wildcards: expand("final_output/{asmname}.full.fa.fai", asmname=get_all_accessions_from_asmset(wildcards.asmset)),
    output:
        report("results/{asmset}/3.quality_assessment/10.ntsynt/{asmset}.k{mink}.w{minw}_ribbon-plot.png",
            category="Quality assessment",
            subcategory="Collinearity",
            caption="../../report/ntsynt.rst",
            labels={"type": "ntSynt", "set": "{asmset}", "k": "{mink}", "w": "{minw}"}),
    log:
        "results/logs/3.quality_assessment/visualise_ntsynt/{asmset}.k{mink}.w{minw}.log"
    benchmark:
        "results/benchmarks/3.quality_assessment/visualise_ntsynt/{asmset}.k{mink}.w{minw}.txt"
    params:
        minlen = 10000000 #minimum length for a block
    container:
        "oras://ghcr.io/dirkjanvw/mogaaap/ntsynt-viz.v1.0.0:latest"
    shell:
        """
        (
        blocks=$(realpath {input.blocks})
        fais="$(realpath {input.fais})"
        cd $(dirname {output})
        ntsynt_viz.py --blocks ${{blocks}} --fais ${{fais}} --seq_length {params.minlen} --prefix $(basename {output} | rev | cut -d '_' -f 2- | rev)
        ) 2> {log}
        """
