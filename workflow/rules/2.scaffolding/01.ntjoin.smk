rule ntjoin:
    input:
        reference = lambda wildcards: config["ref_genome"][wildcards.reference],
        contigs = "results/{asmname}/1.assembly/02.contigs/{asmname}.min{minlen}.sorted.renamed.fa",
    output:
        contigs = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}",
        contigsfai = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.fai",
        contigscount = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{k}.tsv",
        agp = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.agp",
        mxdot = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.mx.dot",
        path = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.path",
        all = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.fa",
        assigned = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.assigned.scaffolds.fa",
        unassigned = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.unassigned.scaffolds.fa",
        unassignedbed = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.tsv.unassigned.bed",
        reference = "results/{asmname}/2.scaffolding/01.ntjoin/{reference}.fa",
        referencefai = "results/{asmname}/2.scaffolding/01.ntjoin/{reference}.fa.fai",
        referencecount = "results/{asmname}/2.scaffolding/01.ntjoin/{reference}.fa.k{k}.w{w}.tsv",
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
        ln -sf $(realpath {input.contigs}) {output.contigs}
        ln -sf $(realpath {input.reference}) {output.reference}
        cd $(dirname {output.all})
        ntJoin assemble target=$(basename {output.contigs}) references='{wildcards.reference}.fa' target_weight='1' reference_weights='2' n=2 G=10000 agp=True no_cut=True overlap=False k={wildcards.k} w={wildcards.w} mkt=True prefix=$(basename {output.all} | rev | cut -d '.' -f 2- | rev) t={threads}
        ) &> {log}
        """
