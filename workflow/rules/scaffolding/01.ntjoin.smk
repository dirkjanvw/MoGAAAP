rule ntjoin:
    input:
        reference = lambda wildcards: config["ref_genome"][wildcards.ref_gen],
        contigs = "results/{asmname}/1.assembly/01.hifiasm/{asmname}.min{minlen}.sorted.renamed.fa"
    output:
        all = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{ref_gen}.min{minlen}.k{k}.w{w}.n1.all.scaffolds.fa",
        assigned = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{ref_gen}.min{minlen}.k{k}.w{w}.n1.assigned.scaffolds.fa",
        unassigned = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{ref_gen}.min{minlen}.k{k}.w{w}.n1.unassigned.scaffolds.fa"
    log:
        "results/logs/ntjoin/{asmname}.vs.{ref_gen}.min{minlen}.k{k}.w{w}.log"
    benchmark:
        "results/benchmarks/ntjoin/{asmname}.vs.{ref_gen}.min{minlen}.k{k}.w{w}.txt"
    threads:
        5
    conda:
        "../../envs/ntjoin.yaml"
    shell:
        """
        (
        ln -s $(realpath {input.contigs}) $(dirname {output.all})/Scaffolds_HiFiasm_{wildcards.asmname}.{wildcards.minlen}
        cd $(dirname {output.all})
        ntJoin assemble target=Scaffolds_HiFiasm_{wildcards.asmname}.vs.{wildcards.ref_gen}.{wildcards.minlen} references='{input.reference}' target_weight='1' reference_weights='2' G=10000 agp=True no_cut=True overlap=False k={wildcards.k} w={wildcards.w} mkt=True prefix=$(basename {output.all} | rev | cut -d '.' -f 2- | rev) t={threads}
        ) &> {log}
        """