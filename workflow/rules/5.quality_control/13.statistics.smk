rule individual_statistics:
    input:
        assembly = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa",
        contigs = expand("results/{{asmname}}/1.assembly/02.contigs/{{asmname}}.min{minlen}.sorted.renamed.fa", minlen=config["min_contig_len"]),
        assigned_sequences = expand("results/{{asmname}}/2.scaffolding/01.ntjoin/{{asmname}}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.assigned.scaffolds.fa", reference=config["ref_genome"], minlen=config["min_contig_len"], k=config["ntjoin_k"], w=config["ntjoin_w"]),
        unassigned_sequences = expand("results/{{asmname}}/2.scaffolding/01.ntjoin/{{asmname}}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.unassigned.scaffolds.fa", reference=config["ref_genome"], minlen=config["min_contig_len"], k=config["ntjoin_k"], w=config["ntjoin_w"]),
        qv = lambda wildcards: expand("results/{{asmname}}/5.quality_control/01.merqury/{k}/{sample}/{{asmname}}_vs_{sample}.{{asmname}}.qv", k=config["k_qc"], sample=config["reads"]["hifi"][wildcards.asmname]),
        reads = lambda wildcards: [config["reads"]["hifi"][wildcards.asmname][sample] for sample in config["reads"]["hifi"][wildcards.asmname]],
        annotation = "results/{asmname}/4.annotation/03.combined/{asmname}.gff",
        query_counts = expand("results/{{asmname}}/3.analysis/06.bcovblp/{query_name}.vs.{{asmname}}.items.circos", query_name=config["prot_queries"]),  ### CIRCOS ITEMS ###
        query_fractions = expand("results/{{asmname}}/3.analysis/07.bcovbln/{query_name}.vs.{{asmname}}.fract.circos", query_name=config["nucl_queries"]),  ### CIRCOS ITEMS ###
    output:
        assembly = "results/{asmname}/5.quality_control/13.statistics/{asmname}.assembly.tsv",
        contigs = "results/{asmname}/5.quality_control/13.statistics/{asmname}.contigs.tsv",
        assigned_sequences = "results/{asmname}/5.quality_control/13.statistics/{asmname}.assigned_sequences.tsv",
        unassigned_sequences = "results/{asmname}/5.quality_control/13.statistics/{asmname}.unassigned_sequences.tsv",
        reads = "results/{asmname}/5.quality_control/13.statistics/{asmname}.reads.tsv",
        tsv = "results/{asmname}/5.quality_control/13.statistics/{asmname}.tsv"
    log:
        "results/logs/5.quality_control/individual_statistics/{asmname}.log"
    benchmark:
        "results/benchmarks/5.quality_control/individual_statistics/{asmname}.txt"
    conda:
        "../../envs/seqkit.yaml"
    shell:
        """
        (
        seqkit stats -abTj1 {input.assembly} > {output.assembly}
        seqkit stats -abTj1 {input.contigs} > {output.contigs}
        seqkit stats -abTj1 {input.assigned_sequences} > {output.assigned_sequences}
        seqkit stats -abTj1 {input.unassigned_sequences} > {output.unassigned_sequences}
        seqkit stats -abTj1 {input.reads} > {output.reads}
        printf "Name\\tTotal length\\t#sequences\\tOverall N50\\t#contigs\\tContig N50\\t#assigned sequences\\t#unassigned sequences\\tTotal length (assigned sequences)\\tTotal length (unassigned sequences)\\tTotal QV (HiF)\\t#HiFi reads\\tN50 HiFi reads\\t#genes\\t#transcripts\\t#queries\\n" > {output.tsv}
        printf "{wildcards.asmname}\\t" >> {output.tsv}
        awk 'BEGIN{{FS = "\\t";}} NR==2{{printf "%s\\t%s\\t%s\\t", $5,$4,$13;}}' {output.assembly} >> {output.tsv}
        awk 'BEGIN{{FS = "\\t";}} NR==2{{printf "%s\\t%s\\t", $4,$13;}}' {output.contigs} >> {output.tsv}
        awk 'BEGIN{{FS = "\\t";}} NR==2{{printf "%s\\t%s\\t", $4,$5;}}' {output.assigned_sequences} >> {output.tsv}
        awk 'BEGIN{{FS = "\\t";}} NR==2{{printf "%s\\t%s\\t", $4,$5;}}' {output.unassigned_sequences} >> {output.tsv}
        awk 'BEGIN{{FS = "\\t";}} NR==1{{printf "%s\\t", $4;}}' {input.qv} >> {output.tsv}
        awk 'BEGIN{{FS = "\\t";}} NR==2{{printf "%s\\t%s\\t", $4,$13;}}' {output.reads} >> {output.tsv}
        awk 'BEGIN{{FS = "\\t";}} $3=="gene"{{genes+=1;}} END{{printf "%s\\t", genes;}}' {input.annotation} >> {output.tsv}
        awk 'BEGIN{{FS = "\\t";}} $3=="mRNA"{{transcripts+=1;}} END{{printf "%s\\t", transcripts;}}' {input.annotation} >> {output.tsv}
        ls {input.query_counts} {input.query_fractions} | wc -l >> {output.tsv}
        ) &> {log}
        """

rule overall_statistics:
    input:
        lambda wildcards: expand("results/{asmname}/5.quality_control/13.statistics/{asmname}.tsv", asmname=config["set"][wildcards.asmset]),
    output:
        "results/{asmset}/5.quality_control/13.statistics/{asmset}.tsv",
    log:
        "results/logs/5.quality_control/overall_statistics/{asmset}.log"
    benchmark:
        "results/benchmarks/5.quality_control/overall_statistics/{asmset}.txt"
    shell:
        "awk 'NR==1 || FNR>1' {input} > {output} 2> {log}"
