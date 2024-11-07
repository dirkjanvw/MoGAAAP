rule individual_statistics:
    input:
        assembly = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa",
        contigs = expand("results/{{asmname}}/1.assembly/02.contigs/{{asmname}}.min{minlen}.sorted.renamed.fa", minlen=config["min_contig_len"]),
        qv = lambda wildcards: expand("results/{{asmname}}/5.quality_assessment/01.merqury/{k}/hifi/{{asmname}}_vs_hifi.qv", k=config["k_qa"]),
        full_annotation = "results/{asmname}/4.annotation/03.combined/{asmname}.gff",
        coding_annotation = "results/{asmname}/4.annotation/03.combined/{asmname}.coding.gff",
    output:
        assembly = "results/{asmname}/5.quality_assessment/13.statistics/{asmname}.assembly.tsv",
        contigs = "results/{asmname}/5.quality_assessment/13.statistics/{asmname}.contigs.tsv",
        chromosomes = "results/{asmname}/5.quality_assessment/13.statistics/{asmname}.assigned_sequences.tsv",
        unassigned = "results/{asmname}/5.quality_assessment/13.statistics/{asmname}.unassigned_sequences.tsv",
        tsv = "results/{asmname}/5.quality_assessment/13.statistics/{asmname}.tsv"
    log:
        "results/logs/5.quality_assessment/individual_statistics/{asmname}.log"
    benchmark:
        "results/benchmarks/5.quality_assessment/individual_statistics/{asmname}.txt"
    params:
        inputdata = lambda wildcards: "HiFi+ONT+Hi-C" if has_ont(wildcards.asmname) and has_hic(wildcards.asmname) else "HiFi+ONT" if has_ont(wildcards.asmname) else "HiFi+Hi-C" if has_hic(wildcards.asmname) else "HiFi only",
        assembler = config["assembler"],
    conda:
        "../../envs/seqkit.yaml"
    shell:
        """
        (
        printf "Name\\tTotal length\\t#sequences\\tN50\\t#genes (full)\\t#genes (coding)\\t#transcripts (coding)\\t#chromosomes\\tTotal length (chromosomes)\\t#unassigned sequences\\tTotal length (unassigned sequences)\\tTotal QV (HiFi)\\t#contigs\\tContig N50\\tInput data\\tAssembler\\n" > {output.tsv}
        printf "{wildcards.asmname}\\t" >> {output.tsv}
        seqkit stats -abTj1 {input.assembly} > {output.assembly}
        awk 'BEGIN{{FS = "\\t";}} NR==2{{printf "%s\\t%s\\t%s\\t", $5,$4,$13;}}' {output.assembly} >> {output.tsv}
        awk 'BEGIN{{FS = "\\t";}} $3=="gene"{{genes+=1;}} END{{printf "%s\\t", genes;}}' {input.full_annotation} >> {output.tsv}
        awk 'BEGIN{{FS = "\\t";}} $3=="gene"{{genes+=1;}} END{{printf "%s\\t", genes;}}' {input.coding_annotation} >> {output.tsv}
        awk 'BEGIN{{FS = "\\t";}} $3=="mRNA"{{transcripts+=1;}} END{{printf "%s\\t", transcripts;}}' {input.coding_annotation} >> {output.tsv}
        seqkit grep -rp "Chr" {input.assembly} | seqkit stats -abTj1 > {output.chromosomes}
        awk 'BEGIN{{FS = "\\t";}} NR==2{{printf "%s\\t%s\\t", $4,$5;}}' {output.chromosomes} >> {output.tsv}
        seqkit grep -vrp "Chr" {input.assembly} | seqkit stats -abTj1 > {output.unassigned}
        awk 'BEGIN{{FS = "\\t";}} NR==2{{printf "%s\\t%s\\t", $4,$5;}}' {output.unassigned} >> {output.tsv}
        awk 'BEGIN{{FS = "\\t";}} NR==1{{printf "%s\\t", $4;}}' {input.qv} >> {output.tsv}
        seqkit stats -abTj1 {input.contigs} > {output.contigs}
        awk 'BEGIN{{FS = "\\t";}} NR==2{{printf "%s\\t%s\\t", $4,$13;}}' {output.contigs} >> {output.tsv}
        printf "{params.inputdata}\\t" >> {output.tsv}
        echo "{params.assembler}" >> {output.tsv}
        ) &> {log}
        """

rule overall_statistics:
    input:
        lambda wildcards: expand("results/{asmname}/5.quality_assessment/13.statistics/{asmname}.tsv", asmname=get_all_accessions_from_asmset(wildcards.asmset, 1)),
    output:
        "results/{asmset}/5.quality_assessment/13.statistics/{asmset}.tsv",
    log:
        "results/logs/5.quality_assessment/overall_statistics/{asmset}.log"
    benchmark:
        "results/benchmarks/5.quality_assessment/overall_statistics/{asmset}.txt"
    shell:
        "awk 'NR==1 || FNR>1' {input} > {output} 2> {log}"

rule visualise_overall_statistics:
    input:
        "results/{asmset}/5.quality_assessment/13.statistics/{asmset}.tsv",
    output:
        report("results/{asmset}/5.quality_assessment/13.statistics/{asmset}.html",
            category="General statistics",
            caption="../../report/statistics.rst",
            labels={"set": "{asmset}"}),
    log:
        "results/logs/5.quality_assessment/visualise_overall_statistics/{asmset}.log"
    benchmark:
        "results/benchmarks/5.quality_assessment/visualise_overall_statistics/{asmset}.txt"
    conda:
        "../../envs/csvtotable.yaml"
    shell:
        "csvtotable -d $'\\t' {input} {output} &> {log}"
