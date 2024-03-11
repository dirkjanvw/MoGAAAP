rule create_renaming_table:
    input:
        agp= "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.all.scaffolds.agp",
        mxdot= "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.all.scaffolds.mx.dot",
    output:
        "results/{asmname}/2.scaffolding/03.renaming/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.all.scaffolds.conversion.tsv"
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

rule renaming_sequences:
    input:
        assigned = expand("results/{{asmname}}/2.scaffolding/01.ntjoin/{{asmname}}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.assigned.scaffolds.fa",
            reference=config["ref_genome"],
            minlen=config["min_contig_len"],
            k=config["ntjoin_k"],
            w=config["ntjoin_w"],
        ),
        unassigned = expand("results/{{asmname}}/2.scaffolding/01.ntjoin/{{asmname}}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.unassigned.scaffolds.fa",
            reference=config["ref_genome"],
            minlen=config["min_contig_len"],
            k=config["ntjoin_k"],
            w=config["ntjoin_w"],
        ),
        table = expand("results/{{asmname}}/2.scaffolding/03.renaming/{{asmname}}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.all.scaffolds.conversion.tsv",
            reference=config["ref_genome"],
            minlen=config["min_contig_len"],
            k=config["ntjoin_k"],
            w=config["ntjoin_w"],
        ),
    output:
        "results/{asmname}/2.scaffolding/03.renaming/{asmname}.fa"
    log:
        "results/logs/2.scaffolding/renaming_sequences/{asmname}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/renaming_sequences/{asmname}.txt"
    params:
        chroms = config["ref_chr"],
    script:
        "../../scripts/renaming_sequences.py"

rule index_sequences:
    input:
        "results/{asmname}/2.scaffolding/03.renaming/{asmname}.fa"
    output:
        "results/{asmname}/2.scaffolding/03.renaming/{asmname}.fa.fai"
    log:
        "results/logs/2.scaffolding/index_sequences/{asmname}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/index_sequences/{asmname}.txt"
    conda:
        "../../envs/samtools.yaml"
    shell:
        "samtools faidx {input} &> {log}"
