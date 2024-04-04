rule ntsynt:
    input:
        genomes = lambda wildcards: expand("results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa", asmname=config["set"][wildcards.asmset]),
        divergence = "results/{asmset}/5.quality_control/09.mash/{asmset}.tsv",
    output:
        blocks = "results/{asmset}/5.quality_control/10.ntsynt/{asmset}.k{mink}.w{minw}.synteny_blocks.tsv",
        commonbf = "results/{asmset}/5.quality_control/10.ntsynt/{asmset}.k{mink}.w{minw}.common.bf",
        mxdot = "results/{asmset}/5.quality_control/10.ntsynt/{asmset}.k{mink}.w{minw}.mx.dot",
        precolblocks = "results/{asmset}/5.quality_control/10.ntsynt/{asmset}.k{mink}.w{minw}.pre-collinear-merge.synteny_blocks.tsv",
    log:
        "results/logs/5.quality_control/ntsynt/{asmset}.k{mink}.w{minw}.log"
    benchmark:
        "results/benchmarks/5.quality_control/ntsynt/{asmset}.k{mink}.w{minw}.txt"
    threads:
        min(workflow.cores, 50)
    conda:
        "../../envs/ntsynt.yaml"
    shell:
        """
        (
        divergence=$(awk '{{for (i=2;i<=NF;i++){{if ($i>d){{d=$i;}}}}}} END{{printf "%.2f",d*100;}}' {input.divergence})
        ln -s $(realpath {input.genomes}) $(dirname {output.blocks})/
        cd $(dirname {output.blocks})
        ntSynt -p {wildcards.asmset}.k{wildcards.mink}.w{wildcards.minw} -k {wildcards.mink} -w {wildcards.minw} -d ${{divergence}} -f -t {threads} $(for file in {input.genomes}; do basename ${{file}}; done)
        ) &> {log}
        """

rule format_ntsynt:
    input:
        fai = lambda wildcards: expand("results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa.fai", asmname=config["set"][wildcards.asmset]),
        blocks = "results/{asmset}/5.quality_control/10.ntsynt/{asmset}.k{mink}.w{minw}.synteny_blocks.tsv",
    output:
        links = "results/{asmset}/5.quality_control/10.ntsynt/{asmset}.k{mink}.w{minw}.links.tsv",
        sequence_lengths = "results/{asmset}/5.quality_control/10.ntsynt/{asmset}.k{mink}.w{minw}.sequence_lengths.tsv",
    log:
        "results/logs/5.quality_control/format_ntsynt/{asmset}.k{mink}.w{minw}.log"
    benchmark:
        "results/benchmarks/5.quality_control/format_ntsynt/{asmset}.k{mink}.w{minw}.txt"
    params:
        minlen = 100000 #minimum length for a block
    conda:
        "../../envs/ntsynt.yaml"
    shell:
        "workflow/scripts/ntSynt.v1.0.0/format_blocks_gggenomes.py -l {params.minlen} -p $(echo {output.links} | rev | cut -d '.' -f 3- | rev) --blocks {input.blocks} --fai {input.fai} &> {log}"

rule visualise_ntsynt:
    input:
        links = "results/{asmset}/5.quality_control/10.ntsynt/{asmset}.k{mink}.w{minw}.links.tsv",
        sequence_lengths = "results/{asmset}/5.quality_control/10.ntsynt/{asmset}.k{mink}.w{minw}.sequence_lengths.tsv",
    output:
        report("results/{asmset}/5.quality_control/10.ntsynt/{asmset}.k{mink}.w{minw}.png",
            category="Collinearity",
            caption="../../report/ntsynt.rst",
            labels={"type": "ntSynt", "set": "{asmset}", "k": "{mink}", "w": "{minw}"}),
    log:
        "results/logs/5.quality_control/visualise_ntsynt/{asmset}.k{mink}.w{minw}.log"
    benchmark:
        "results/benchmarks/5.quality_control/visualise_ntsynt/{asmset}.k{mink}.w{minw}.txt"
    params:
        minlen = 10000000 #minimum length for a block
    container:
        "workflow/singularity/ntSynt/ntSynt.visualization_scripts.v1.0.0.sif"
    shell:
        "workflow/scripts/ntSynt.v1.0.0/plot_synteny_blocks_gggenomes.R -s {input.sequence_lengths} -l {input.links} -p $(echo {output} | rev | cut -d '.' -f 2- | rev) &> {log}"