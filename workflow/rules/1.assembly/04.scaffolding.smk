rule ntjoin:
    input:
        reference = lambda wildcards: config["reference_genomes"][wildcards.reference]["genome"],
        contigs = "results/{asmname}/1.assembly/02.contigs/{asmname}.min{minlen}.sorted.renamed.fa",
    output:
        all = "results/{asmname}/1.assembly/04.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.fa",
        assigned = "results/{asmname}/1.assembly/04.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.assigned.scaffolds.fa",
        unassigned = "results/{asmname}/1.assembly/04.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.unassigned.scaffolds.fa",
        agp = "results/{asmname}/1.assembly/04.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.agp",
        mxdot = "results/{asmname}/1.assembly/04.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.mx.dot",
        contigscount = "results/{asmname}/1.assembly/04.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.tsv",
        path = "results/{asmname}/1.assembly/04.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.path",
        unassignedbed = "results/{asmname}/1.assembly/04.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.tsv.unassigned.bed",
    log:
        "results/logs/1.assembly/ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.log"
    benchmark:
        "results/benchmarks/1.assembly/ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.txt"
    threads:
        5
    conda:
        "../../envs/ntjoin.yaml"
    shell:
        """
        (
        ln -sf $(realpath {input.contigs}) $(dirname {output.all})/$(basename {output.all} | rev | cut -d"." -f7- | rev)
        ln -sf $(realpath {input.reference}) $(dirname {output.all})/{wildcards.reference}.fa
        cd $(dirname {output.all})
        ntJoin assemble target=$(basename {output.all} | rev | cut -d"." -f7- | rev) references='{wildcards.reference}.fa' target_weight='1' reference_weights='2' n=2 G=10000 agp=True no_cut=True overlap=False k={wildcards.k} w={wildcards.w} mkt=True prefix=$(basename {output.all} | rev | cut -d '.' -f 2- | rev) t={threads}
        ) &> {log}
        """

rule ragtag:
    input:
        reference = lambda wildcards: config["reference_genomes"][wildcards.reference]["genome"],
        contigs = "results/{asmname}/1.assembly/02.contigs/{asmname}.min{minlen}.sorted.renamed.fa",
    output:
        agp = "results/{asmname}/1.assembly/04.ragtag/{asmname}.vs.{reference}.min{minlen}/ragtag.scaffold.agp",
        paf = "results/{asmname}/1.assembly/04.ragtag/{asmname}.vs.{reference}.min{minlen}/ragtag.scaffold.asm.paf",
        paflog = "results/{asmname}/1.assembly/04.ragtag/{asmname}.vs.{reference}.min{minlen}/ragtag.scaffold.asm.paf.log",
        confidence = "results/{asmname}/1.assembly/04.ragtag/{asmname}.vs.{reference}.min{minlen}/ragtag.scaffold.confidence.txt",
        err = "results/{asmname}/1.assembly/04.ragtag/{asmname}.vs.{reference}.min{minlen}/ragtag.scaffold.err",
        fasta = "results/{asmname}/1.assembly/04.ragtag/{asmname}.vs.{reference}.min{minlen}/ragtag.scaffold.fasta",
        stats = "results/{asmname}/1.assembly/04.ragtag/{asmname}.vs.{reference}.min{minlen}/ragtag.scaffold.stats",
    log:
        "results/logs/1.assembly/ragtag/{asmname}.vs.{reference}.min{minlen}.log"
    benchmark:
        "results/benchmarks/1.assembly/ragtag/{asmname}.vs.{reference}.min{minlen}.txt"
    threads:
        5
    conda:
        "../../envs/ragtag.yaml"
    shell:
        "ragtag.py scaffold -rt{threads} -o$(dirname {output.fasta}) {input.reference} {input.contigs} &> {log}"

rule bwa_index_contigs:
    input:
        contigs = "results/{asmname}/1.assembly/02.contigs/{asmname}.min{minlen}.sorted.renamed.fa",
    output:
        index1 = "results/{asmname}/1.assembly/02.contigs/{asmname}.min{minlen}.sorted.renamed.fa.0123",
        index2 = "results/{asmname}/1.assembly/02.contigs/{asmname}.min{minlen}.sorted.renamed.fa.amb",
        index3 = "results/{asmname}/1.assembly/02.contigs/{asmname}.min{minlen}.sorted.renamed.fa.ann",
        index4 = "results/{asmname}/1.assembly/02.contigs/{asmname}.min{minlen}.sorted.renamed.fa.bwt.2bit.64",
        index5 = "results/{asmname}/1.assembly/02.contigs/{asmname}.min{minlen}.sorted.renamed.fa.pac",
    log:
        "results/logs/1.assembly/bwa_index_contigs/{asmname}.min{minlen}.log"
    benchmark:
        "results/benchmarks/1.assembly/bwa_index_contigs/{asmname}.min{minlen}.txt"
    conda:
        "../../envs/mapping.yaml"
    shell:
        "bwa-mem2 index {input.contigs} &> {log}"

rule map_hic:
    input:
        contigs = expand("results/{{asmname}}/1.assembly/02.contigs/{{asmname}}.min{minlen}.sorted.renamed.fa", minlen=config["min_contig_len"],),
        index1 = expand("results/{{asmname}}/1.assembly/02.contigs/{{asmname}}.min{minlen}.sorted.renamed.fa.0123", minlen=config["min_contig_len"],),
        index2 = expand("results/{{asmname}}/1.assembly/02.contigs/{{asmname}}.min{minlen}.sorted.renamed.fa.amb", minlen=config["min_contig_len"],),
        index3 = expand("results/{{asmname}}/1.assembly/02.contigs/{{asmname}}.min{minlen}.sorted.renamed.fa.ann", minlen=config["min_contig_len"],),
        index4 = expand("results/{{asmname}}/1.assembly/02.contigs/{{asmname}}.min{minlen}.sorted.renamed.fa.bwt.2bit.64", minlen=config["min_contig_len"],),
        index5 = expand("results/{{asmname}}/1.assembly/02.contigs/{{asmname}}.min{minlen}.sorted.renamed.fa.pac", minlen=config["min_contig_len"],),
        forward = get_hic_1,
        backward = get_hic_2,
    output:
        "results/{asmname}/1.assembly/04.hic/{asmname}.hic.sorted.bam",
    log:
        "results/logs/1.assembly/map_hic/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/map_hic/{asmname}.txt"
    threads:
        10
    conda:
        "../../envs/mapping.yaml"
    shell:
        "(bwa-mem2 mem -t {threads} -5SP {input.contigs} {input.forward} {input.backward} | samblaster | samtools view - -@ $(({threads}-1)) -S -h -b -F 3340 -o {output}) &> {log}"

rule filter_hic:
    input:
        "results/{asmname}/1.assembly/04.hic/{asmname}.hic.sorted.bam",
    output:
        "results/{asmname}/1.assembly/04.hic/{asmname}.hic.sorted.filtered.bam",
    log:
        "results/logs/1.assembly/filter_hic/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/filter_hic/{asmname}.txt"
    threads:
        10
    container:
        "oras://ghcr.io/dirkjanvw/mogaaap/haphic.f8f7451:latest"
    shell:
        "(filter_bam.py {input} 1 --NM 3 --threads {threads} | samtools view - -b -@ $(({threads}-1)) -o {output}) &> {log}"

rule ntjoin_plot_hic:
    input:
        agp = lambda wildcards: expand("results/{{asmname}}/1.assembly/04.ntjoin/{{asmname}}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.agp",
            reference=get_reference_id(wildcards.asmname),
            minlen=config["min_contig_len"],
            k=config["ntjoin_k"],
            w=config["ntjoin_w"],
        ),
        bam = "results/{asmname}/1.assembly/04.hic/{asmname}.hic.sorted.filtered.bam",
        table = "results/{asmname}/1.assembly/05.renaming/{asmname}.ntjoin.html", #making sure that the table is generated before the plot so that the report can be interpreted
    output:
        pdf = report("results/{asmname}/1.assembly/04.ntjoin/contact_map.pdf",
            category="assembly",
            subcategory="Hi-C",
            caption="../../report/hic.rst",
            labels={"assembly": "{asmname}",
                    "stage": "scaffolds",
                    "algorithm": "ntjoin (contact map)"}
        ),
        pkl = "results/{asmname}/1.assembly/04.ntjoin/contact_matrix.pkl",
    log:
        "results/logs/1.assembly/ntjoin_plot/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/ntjoin_plot/{asmname}.txt"
    threads:
        8
    container:
        "oras://ghcr.io/dirkjanvw/mogaaap/haphic.f8f7451:latest"
    shell:
        """
        (
        agp="$(realpath {input.agp})"
        bam="$(realpath {input.bam})"
        cd $(dirname {output.pdf})
        haphic plot --threads {threads} $agp $bam
        ) &> {log}
        """

use rule ntjoin_plot_hic as ragtag_plot_hic with:
    input:
        agp = lambda wildcards: expand("results/{{asmname}}/1.assembly/04.ragtag/{{asmname}}.vs.{reference}.min{minlen}/ragtag.scaffold.agp",
            reference=get_reference_id(wildcards.asmname),
            minlen=config["min_contig_len"],
        ),
        bam = "results/{asmname}/1.assembly/04.hic/{asmname}.hic.sorted.filtered.bam",
        table = "results/{asmname}/1.assembly/05.renaming/{asmname}.ragtag.html", #making sure that the table is generated before the plot so that the report can be interpreted
    output:
        pdf = report("results/{asmname}/1.assembly/04.ragtag/contact_map.pdf",
            category="assembly",
            subcategory="Hi-C",
            caption="../../report/hic.rst",
            labels={"assembly": "{asmname}",
                    "stage": "scaffolds",
                    "algorithm": "ragtag (contact map)"}
        ),
        pkl = "results/{asmname}/1.assembly/04.ragtag/contact_matrix.pkl",
    log:
        "results/logs/1.assembly/ragtag_plot/{asmname}.log"
    benchmark:
        "results/benchmarks/1.assembly/ragtag_plot/{asmname}.txt"
