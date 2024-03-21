rule create_renaming_table:
    input:
        agp= "results/{asmname}/2.scaffolding/01.ntjoin/k{k}/w{w}/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.all.scaffolds.agp",
        mxdot= "results/{asmname}/2.scaffolding/01.ntjoin/k{k}/w{w}/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.all.scaffolds.mx.dot",
    output:
        "results/{asmname}/2.scaffolding/02.renaming/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.all.scaffolds.conversion.tsv"
    log:
        "results/logs/2.scaffolding/create_renaming_table/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.all.scaffolds.log"
    benchmark:
        "results/benchmarks/2.scaffolding/create_renaming_table/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.all.scaffolds.txt"
    params:
        num_chr = len(config["ref_chr"]),
    shell: # I found the following awk commands to get a quick list of what chromosomes belong together (based on https://github.com/bcgsc/ntJoin/issues/63):
        """
        (
        awk 'BEGIN{{OFS = "'\\''";}} {{print $1,$6;}}' {input.agp} | \
            awk -vasm="{wildcards.asmname}" -vref="{wildcards.reference}" 'BEGIN{{FS = "'\\''"; OFS = "\\t";}} FNR==NR{{chr[$2]=$1; next;}} $1~"^"ref{{printf "%s\\t",$2;}} $1~"^"asm{{printf "%s\\n",chr[$2];}}' - {input.mxdot} | \
            sort | \
            uniq -c | \
            sort -nr | \
            awk -vnum_chr="{params.num_chr}" 'BEGIN{{OFS = "\\t";}} FNR<=num_chr{{print $2,$3;}}' > {output}
        ) &> {log}
        """

rule renaming_scaffolds:
    input:
        assigned = expand("results/{{asmname}}/2.scaffolding/01.ntjoin/k{k}/w{w}/{{asmname}}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.assigned.scaffolds.fa",
            reference=config["ref_genome"],
            minlen=config["min_contig_len"],
            k=config["ntjoin_k"],
            w=config["ntjoin_w"],
        ),
        unassigned = expand("results/{{asmname}}/2.scaffolding/01.ntjoin/k{k}/w{w}/{{asmname}}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.unassigned.scaffolds.fa",
            reference=config["ref_genome"],
            minlen=config["min_contig_len"],
            k=config["ntjoin_k"],
            w=config["ntjoin_w"],
        ),
        table = expand("results/{{asmname}}/2.scaffolding/02.renaming/{{asmname}}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.all.scaffolds.conversion.tsv",
            reference=config["ref_genome"],
            minlen=config["min_contig_len"],
            k=config["ntjoin_k"],
            w=config["ntjoin_w"],
        ),
    output:
        "results/{asmname}/2.scaffolding/02.renaming/{asmname}.unsorted.fa"
    log:
        "results/logs/2.scaffolding/renaming_scaffolds/{asmname}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/renaming_scaffolds/{asmname}.txt"
    params:
        chroms = config["ref_chr"],
    script:
        "../../scripts/renaming_scaffolds.py"

rule sort_scaffolds:
    input:
        "results/{asmname}/2.scaffolding/02.renaming/{asmname}.unsorted.fa"
    output:
        protected("results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa")
    log:
        "results/logs/2.scaffolding/sort_scaffolds/{asmname}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/sort_scaffolds/{asmname}.txt"
    conda:
        "../../envs/seqkit.yaml"
    shell:
        "seqkit sort -n {input} > {output} 2> {log}"

rule index_scaffolds:
    input:
        "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa"
    output:
        "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa.fai"
    log:
        "results/logs/2.scaffolding/index_scaffolds/{asmname}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/index_scaffolds/{asmname}.txt"
    conda:
        "../../envs/samtools.yaml"
    shell:
        "samtools faidx {input} &> {log}"
