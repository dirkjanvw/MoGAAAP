rule individual_statistics_full:
    input:
        assembly = "final_output/{asmname}.full.fa",
        contigs = expand("results/{{asmname}}/1.assembly/02.contigs/{{asmname}}.min{minlen}.sorted.renamed.fa", minlen=config["min_contig_len"]),
        qv = lambda wildcards: expand("results/{{asmname}}/3.quality_assessment/01.merqury/{k}/{wgstype}/{{asmname}}_vs_{wgstype}.qv", k=config["k_qa"], wgstype=get_best_wgstype(wildcards.asmname).lower()),
        kmerstats = lambda wildcards: expand("results/{{asmname}}/3.quality_assessment/01.merqury/{k}/{wgstype}/{{asmname}}_vs_{wgstype}.completeness.stats", k=config["k_qa"], wgstype=get_best_wgstype(wildcards.asmname).lower()),
        full_annotation = "final_output/{asmname}.full.gff",
        coding_annotation = "final_output/{asmname}.full.coding.gff",
    output:
        assembly = "results/{asmname}/3.quality_assessment/13.statistics/{asmname}.assembly.tsv",
        contigs = "results/{asmname}/3.quality_assessment/13.statistics/{asmname}.contigs.tsv",
        chromosomes = "results/{asmname}/3.quality_assessment/13.statistics/{asmname}.assigned_sequences.tsv",
        unassigned = "results/{asmname}/3.quality_assessment/13.statistics/{asmname}.unassigned_sequences.tsv",
        tsv = "results/{asmname}/3.quality_assessment/13.statistics/{asmname}.full.tsv"
    log:
        "results/logs/3.quality_assessment/individual_statistics_full/{asmname}.log"
    benchmark:
        "results/benchmarks/3.quality_assessment/individual_statistics_full/{asmname}.txt"
    params:
        inputdata = lambda wildcards: "HiFi+ONT+Hi-C" if has_ont(wildcards.asmname) and has_hic(wildcards.asmname) else "HiFi+ONT" if has_ont(wildcards.asmname) else "HiFi+Hi-C" if has_hic(wildcards.asmname) else "HiFi only",
        assembler = config["assembler"],
        best_wgstype = lambda wildcards: get_best_wgstype(wildcards.asmname)
    conda:
        "../../envs/seqkit.yaml"
    shell:
        """
        (
        printf "Name\\tTotal length\\t#sequences\\tN50\\tTotal QV\\tK-mer completeness\\t#genes (full)\\t#genes (coding)\\t#transcripts (coding)\\t#chromosomes\\tTotal length (chromosomes)\\t#unassigned sequences\\tTotal length (unassigned sequences)\\t#contigs\\tContig N50\\tInput data\\tAssembler\\n" > {output.tsv}
        printf "{wildcards.asmname}\\t" >> {output.tsv}
        seqkit stats -abTj1 {input.assembly} > {output.assembly}
        awk 'BEGIN{{FS = "\\t";}} NR==2{{printf "%s\\t%s\\t%s\\t", $5,$4,$13;}}' {output.assembly} >> {output.tsv}
        awk 'BEGIN{{FS = "\\t";}} NR==1{{printf "%s ({params.best_wgstype})\\t", $4;}}' {input.qv} >> {output.tsv}
        awk 'BEGIN{{FS = "\\t";}} NR==1{{printf "%s ({params.best_wgstype})\\t", $5;}}' {input.kmerstats} >> {output.tsv}
        awk 'BEGIN{{FS = "\\t";}} $3=="gene"{{genes+=1;}} END{{printf "%s\\t", genes;}}' {input.full_annotation} >> {output.tsv}
        awk 'BEGIN{{FS = "\\t";}} $3=="gene"{{genes+=1;}} END{{printf "%s\\t", genes;}}' {input.coding_annotation} >> {output.tsv}
        awk 'BEGIN{{FS = "\\t";}} $3=="mRNA"{{transcripts+=1;}} END{{printf "%s\\t", transcripts;}}' {input.coding_annotation} >> {output.tsv}
        seqkit grep -rp "Chr" {input.assembly} | seqkit stats -abTj1 > {output.chromosomes}
        awk 'BEGIN{{FS = "\\t";}} NR==2{{printf "%s\\t%s\\t", $4,$5;}}' {output.chromosomes} >> {output.tsv}
        seqkit grep -vrp "Chr" {input.assembly} | seqkit stats -abTj1 > {output.unassigned}
        awk 'BEGIN{{FS = "\\t";}} NR==2{{printf "%s\\t%s\\t", $4,$5;}}' {output.unassigned} >> {output.tsv}
        seqkit stats -abTj1 {input.contigs} > {output.contigs}
        awk 'BEGIN{{FS = "\\t";}} NR==2{{printf "%s\\t%s\\t", $4,$13;}}' {output.contigs} >> {output.tsv}
        printf "{params.inputdata}\\t" >> {output.tsv}
        echo "{params.assembler}" >> {output.tsv}
        ) &> {log}
        """

rule individual_statistics_no_assembly_with_merqury:
    input:
        assembly = "final_output/{asmname}.full.fa",
        qv = lambda wildcards: expand("results/{{asmname}}/3.quality_assessment/01.merqury/{k}/{wgstype}/{{asmname}}_vs_{wgstype}.qv", k=config["k_qa"], wgstype=get_best_wgstype(wildcards.asmname).lower()),
        kmerstats = lambda wildcards: expand("results/{{asmname}}/3.quality_assessment/01.merqury/{k}/{wgstype}/{{asmname}}_vs_{wgstype}.completeness.stats", k=config["k_qa"], wgstype=get_best_wgstype(wildcards.asmname).lower()),
        full_annotation = "final_output/{asmname}.full.gff",
        coding_annotation = "final_output/{asmname}.full.coding.gff",
    output:
        assembly = "results/{asmname}/3.quality_assessment/13.statistics/{asmname}.assembly.tsv",
        tsv = "results/{asmname}/3.quality_assessment/13.statistics/{asmname}.medium.tsv"
    log:
        "results/logs/3.quality_assessment/individual_statistics_no_assembly/{asmname}.log"
    benchmark:
        "results/benchmarks/3.quality_assessment/individual_statistics_no_assembly/{asmname}.txt"
    params:
        best_wgstype = lambda wildcards: get_best_wgstype(wildcards.asmname)
    conda:
        "../../envs/seqkit.yaml"
    shell:
        """
        (
        printf "Name\\tTotal length\\t#sequences\\tN50\\tTotal QV\\tK-mer completeness\\t#genes (full)\\t#genes (coding)\\t#transcripts (coding)\\n" > {output.tsv}
        printf "{wildcards.asmname}\\t" >> {output.tsv}
        seqkit stats -abTj1 {input.assembly} > {output.assembly}
        awk 'BEGIN{{FS = "\\t";}} NR==2{{printf "%s\\t%s\\t%s\\t", $5,$4,$13;}}' {output.assembly} >> {output.tsv}
        awk 'BEGIN{{FS = "\\t";}} NR==1{{printf "%s ({params.best_wgstype})\\t", $4;}}' {input.qv} >> {output.tsv}
        awk 'BEGIN{{FS = "\\t";}} NR==1{{printf "%s ({params.best_wgstype})\\t", $5;}}' {input.kmerstats} >> {output.tsv}
        awk 'BEGIN{{FS = "\\t";}} $3=="gene"{{genes+=1;}} END{{printf "%s\\t", genes;}}' {input.full_annotation} >> {output.tsv}
        awk 'BEGIN{{FS = "\\t";}} $3=="gene"{{genes+=1;}} END{{printf "%s\\t", genes;}}' {input.coding_annotation} >> {output.tsv}
        awk 'BEGIN{{FS = "\\t";}} $3=="mRNA"{{transcripts+=1;}} END{{printf "%s\\t", transcripts;}}' {input.coding_annotation} >> {output.tsv}
        ) &> {log}
        """

rule individual_statistics_no_assembly:
    input:
        assembly = "final_output/{asmname}.full.fa",
        full_annotation = "final_output/{asmname}.full.gff",
        coding_annotation = "final_output/{asmname}.full.coding.gff",
    output:
        assembly = "results/{asmname}/3.quality_assessment/13.statistics/{asmname}.assembly.tsv",
        tsv = "results/{asmname}/3.quality_assessment/13.statistics/{asmname}.small.tsv"
    log:
        "results/logs/3.quality_assessment/individual_statistics_no_assembly/{asmname}.log"
    benchmark:
        "results/benchmarks/3.quality_assessment/individual_statistics_no_assembly/{asmname}.txt"
    conda:
        "../../envs/seqkit.yaml"
    shell:
        """
        (
        printf "Name\\tTotal length\\t#sequences\\tN50\\t#genes (full)\\t#genes (coding)\\t#transcripts (coding)\\n" > {output.tsv}
        printf "{wildcards.asmname}\\t" >> {output.tsv}
        seqkit stats -abTj1 {input.assembly} > {output.assembly}
        awk 'BEGIN{{FS = "\\t";}} NR==2{{printf "%s\\t%s\\t%s\\t", $5,$4,$13;}}' {output.assembly} >> {output.tsv}
        awk 'BEGIN{{FS = "\\t";}} $3=="gene"{{genes+=1;}} END{{printf "%s\\t", genes;}}' {input.full_annotation} >> {output.tsv}
        awk 'BEGIN{{FS = "\\t";}} $3=="gene"{{genes+=1;}} END{{printf "%s\\t", genes;}}' {input.coding_annotation} >> {output.tsv}
        awk 'BEGIN{{FS = "\\t";}} $3=="mRNA"{{transcripts+=1;}} END{{printf "%s\\t", transcripts;}}' {input.coding_annotation} >> {output.tsv}
        ) &> {log}
        """

def get_invididual_statistics(wildcards):
    filelist = []
    for asmname in get_all_accessions_from_asmset(wildcards.asmset, 1):
        if has_assembly_location(asmname):
            if has_illumina(asmname) or has_hifi(asmname):
                filelist.append(f"results/{asmname}/3.quality_assessment/13.statistics/{asmname}.medium.tsv")
            else:
                filelist.append(f"results/{asmname}/3.quality_assessment/13.statistics/{asmname}.small.tsv")
        else:
            filelist.append(f"results/{asmname}/3.quality_assessment/13.statistics/{asmname}.full.tsv")
    return filelist

rule overall_statistics:
    input:
        get_invididual_statistics,
    output:
        "results/{asmset}/3.quality_assessment/13.statistics/{asmset}.tsv",
    log:
        "results/logs/3.quality_assessment/overall_statistics/{asmset}.log"
    benchmark:
        "results/benchmarks/3.quality_assessment/overall_statistics/{asmset}.txt"
    shell:
        "./workflow/scripts/merge_tsv_files.sh {input} > {output} 2> {log}"

rule visualise_overall_statistics:
    input:
        "results/{asmset}/3.quality_assessment/13.statistics/{asmset}.tsv",
    output:
        report("results/{asmset}/3.quality_assessment/13.statistics/{asmset}.html",
            category="Quality assessment",
            subcategory="General statistics",
            caption="../../report/statistics.rst",
            labels={"set": "{asmset}"}),
    log:
        "results/logs/3.quality_assessment/visualise_overall_statistics/{asmset}.log"
    benchmark:
        "results/benchmarks/3.quality_assessment/visualise_overall_statistics/{asmset}.txt"
    conda:
        "../../envs/csvtotable.yaml"
    shell:
        "csvtotable -d $'\\t' {input} {output} &> {log}"
