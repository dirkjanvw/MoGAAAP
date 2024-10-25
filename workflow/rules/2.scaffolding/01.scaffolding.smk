rule ntjoin:
    input:
        reference = lambda wildcards: config["reference_genomes"][wildcards.reference]["genome"],
        contigs = "results/{asmname}/1.assembly/02.contigs/{asmname}.min{minlen}.sorted.renamed.fa",
    output:
        all = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.fa",
        assigned = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.assigned.scaffolds.fa",
        unassigned = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.unassigned.scaffolds.fa",
        agp = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.agp",
        mxdot = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.mx.dot",
        contigscount = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.tsv",
        path = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.path",
        unassignedbed = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.tsv.unassigned.bed",
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
        ln -sf $(realpath {input.contigs}) $(dirname {output.all})/$(basename {output.all} | rev | cut -d"." -f7- | rev)
        ln -sf $(realpath {input.reference}) $(dirname {output.all})/{wildcards.reference}.fa
        cd $(dirname {output.all})
        ntJoin assemble target=$(basename {output.all} | rev | cut -d"." -f7- | rev) references='{wildcards.reference}.fa' target_weight='1' reference_weights='2' n=2 G=10000 agp=True no_cut=True overlap=False k={wildcards.k} w={wildcards.w} mkt=True prefix=$(basename {output.all} | rev | cut -d '.' -f 2- | rev) t={threads}
        ) &> {log}
        """

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
        "results/logs/2.scaffolding/bwa_index_contigs/{asmname}.min{minlen}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/bwa_index_contigs/{asmname}.min{minlen}.txt"
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
        "results/{asmname}/2.scaffolding/01.hic/{asmname}.hic.sorted.bam",
    log:
        "results/logs/2.scaffolding/map_hic/{asmname}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/map_hic/{asmname}.txt"
    threads:
        10
    conda:
        "../../envs/mapping.yaml"
    shell:
        "(bwa-mem2 mem -t {threads} -5SP {input.contigs} {input.forward} {input.backward} | samblaster | samtools view - -@ $(({threads}-1)) -S -h -b -F 3340 -o {output}) &> {log}"

rule filter_hic:
    input:
        "results/{asmname}/2.scaffolding/01.hic/{asmname}.hic.sorted.bam",
    output:
        "results/{asmname}/2.scaffolding/01.hic/{asmname}.hic.sorted.filtered.bam",
    log:
        "results/logs/2.scaffolding/filter_hic/{asmname}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/filter_hic/{asmname}.txt"
    threads:
        10
    singularity:
        "workflow/singularity/haphic/haphic.f8f7451.sif"
    shell:
        "(filter_bam.py {input} 1 --NM 3 --threads {threads} | samtools view - -b -@ $(({threads}-1)) -o {output}) &> {log}"

rule haphic:
    input:
        contigs = expand("results/{{asmname}}/1.assembly/02.contigs/{{asmname}}.min{minlen}.sorted.renamed.fa", minlen=config["min_contig_len"],),
        hic = "results/{asmname}/2.scaffolding/01.hic/{asmname}.hic.sorted.filtered.bam",
        #gfa = "", #TODO: we could at some point look into using GFA from hifiasm if hifiasm was used, but this interferes with the current renaming of contigs (HapHiC calls this feature EXPERIMENTAL!)
    output:
        fa = "results/{asmname}/2.scaffolding/01.haphic/{asmname}_HapHiC/04.build/scaffolds.fa",
        salse_agp = "results/{asmname}/2.scaffolding/01.haphic/{asmname}_HapHiC/04.build/scaffolds.agp",
        yahs_agp = "results/{asmname}/2.scaffolding/01.haphic/{asmname}_HapHiC/04.build/scaffolds.raw.agp",
    log:
        "results/logs/2.scaffolding/haphic/{asmname}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/haphic/{asmname}.txt"
    params:
        num_chr = lambda wildcards: len(config["reference_genomes"][get_reference_id(wildcards.asmname)]["chromosomes"]),  #the number of chromosomes
    threads:
        8
    singularity:
        "workflow/singularity/haphic/haphic.f8f7451.sif"
    shell:
        "haphic pipeline --threads {threads} --outdir $(dirname $(dirname {output.fa})) --verbose {input.contigs} {input.hic} {params.num_chr} &> {log}"

rule haphic_plot:
    input:
        agp = "results/{asmname}/2.scaffolding/01.haphic/{asmname}_HapHiC/04.build/scaffolds.raw.agp",
        bam = "results/{asmname}/2.scaffolding/01.hic/{asmname}.hic.sorted.filtered.bam",
    output:
        pdf = report("results/{asmname}/2.scaffolding/01.haphic/contact_map.pdf",
            category="Hi-C",
            caption="../../report/hic.rst",
            labels={"assembly": "{asmname}",
                    "algorithm": "HapHiC"}
        ),
        pkl = "results/{asmname}/2.scaffolding/01.haphic/contact_matrix.pkl",
    log:
        "results/logs/2.scaffolding/haphic_plot/{asmname}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/haphic_plot/{asmname}.txt"
    threads:
        8
    singularity:
        "workflow/singularity/haphic/haphic.f8f7451.sif"
    shell:
        """
        (
        agp="$(realpath {input.agp})"
        bam="$(realpath {input.bam})"
        cd $(dirname {output.pdf})
        haphic plot --threads {threads} $agp $bam
        ) &> {log}
        """

rule yahs:
    input:
        contigs = expand("results/{{asmname}}/1.assembly/02.contigs/{{asmname}}.min{minlen}.sorted.renamed.fa", minlen=config["min_contig_len"],),
        contigs_fai = expand("results/{{asmname}}/1.assembly/02.contigs/{{asmname}}.min{minlen}.sorted.renamed.fa.fai", minlen=config["min_contig_len"],),
        hic = "results/{asmname}/2.scaffolding/01.hic/{asmname}.hic.sorted.filtered.bam",
    output:
        fa = "results/{asmname}/2.scaffolding/01.yahs/{asmname}_scaffolds_final.fa",
        agp = "results/{asmname}/2.scaffolding/01.yahs/{asmname}_scaffolds_final.agp",
    log:
        "results/logs/2.scaffolding/yahs/{asmname}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/yahs/{asmname}.txt"
    params:
        telomere_motif = config["telomere_motif"],
    threads:
        8
    conda:
        "../../envs/yahs.yaml"
    shell:
        """
        (
        mkdir -p $(dirname {output.fa})
        yahs --telo-motif {params.telomere_motif} -o $(echo {output.fa} | sed 's/_scaffolds_final.fa$//g') {input.contigs} {input.hic}
        ) &> {log}
        """

use rule haphic_plot as yahs_plot with:
    input:
        agp = "results/{asmname}/2.scaffolding/01.yahs/{asmname}_scaffolds_final.agp",
        bam = "results/{asmname}/2.scaffolding/01.hic/{asmname}.hic.sorted.filtered.bam",
    output:
        pdf = report("results/{asmname}/2.scaffolding/01.yahs/contact_map.pdf",
            category="Hi-C",
            caption="../../report/hic.rst",
            labels={"assembly": "{asmname}",
                    "algorithm": "YaHS"}
        ),
        pkl = "results/{asmname}/2.scaffolding/01.yahs/contact_matrix.pkl",
    log:
        "results/logs/2.scaffolding/yahs_plot/{asmname}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/yahs_plot/{asmname}.txt"

use rule haphic_plot as ntjoin_plot with:
    input:
        agp = lambda wildcards: expand("results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.agp",
            reference=get_reference_id(wildcards.asmname),
            minlen=config["min_contig_len"],
            k=config["ntjoin_k"],
            w=config["ntjoin_w"],
        ),
        bam = "results/{asmname}/2.scaffolding/01.hic/{asmname}.hic.sorted.filtered.bam",
    output:
        pdf = report("results/{asmname}/2.scaffolding/01.ntjoin/contact_map.pdf",
            category="Hi-C",
            caption="../../report/hic.rst",
            labels={"assembly": "{asmname}",
                    "algorithm": "ntJoin"}
        ),
        pkl = "results/{asmname}/2.scaffolding/01.ntjoin/contact_matrix.pkl",
    log:
        "results/logs/2.scaffolding/ntjoin_plot/{asmname}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/ntjoin_plot/{asmname}.txt"








































