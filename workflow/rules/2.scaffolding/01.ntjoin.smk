rule ntjoin:
    input:
        reference = "results/{asmname}/2.scaffolding/00.reference/{reference}.fa",
        contigs = "results/{asmname}/1.assembly/01.hifiasm/{asmname}.min{minlen}.sorted.renamed.fa",
    output:
        all = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.all.scaffolds.fa",
        assigned = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.assigned.scaffolds.fa",
        unassigned = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n1.unassigned.scaffolds.fa",
    log:
        "results/logs/2.scaffolding/ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.txt"
    threads:
        5
    conda:
        "../../envs/ntjoin.yaml"
    shell:
        """
        (
        ln -s $(realpath {input}) $(dirname {output.all})/
        cd $(dirname {output.all})
        ntJoin assemble target=$(basename {input.contigs}) references='$(basename {input.reference})' target_weight='1' reference_weights='2' G=10000 agp=True no_cut=True overlap=False k={wildcards.k} w={wildcards.w} mkt=True prefix=$(basename {output.all} | rev | cut -d '.' -f 2- | rev) t={threads}
        ) &> {log}
        """
