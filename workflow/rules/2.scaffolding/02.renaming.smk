rule create_renaming_table:
    input:
        agp= "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.agp",
        mxdot= "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.mx.dot",
    output:
        "results/{asmname}/2.scaffolding/02.renaming/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.conversion.tsv"
    log:
        "results/logs/2.scaffolding/create_renaming_table/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.log"
    benchmark:
        "results/benchmarks/2.scaffolding/create_renaming_table/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.txt"
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
        all = expand("results/{{asmname}}/2.scaffolding/01.ntjoin/{{asmname}}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.fa",
            reference=config["ref_genome"],
            minlen=config["min_contig_len"],
            k=config["ntjoin_k"],
            w=config["ntjoin_w"],
        ),
        table = expand("results/{{asmname}}/2.scaffolding/02.renaming/{{asmname}}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.conversion.tsv",
            reference=config["ref_genome"],
            minlen=config["min_contig_len"],
            k=config["ntjoin_k"],
            w=config["ntjoin_w"],
        ),
    output:
        "results/{asmname}/2.scaffolding/02.renaming/{asmname}.unsorted.unoriented.fa"
    log:
        "results/logs/2.scaffolding/renaming_scaffolds/{asmname}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/renaming_scaffolds/{asmname}.txt"
    params:
        chroms = config["ref_chr"],
    script:
        "../../scripts/renaming_scaffolds.py"

rule obtain_chromosomes:
    input:
        "results/{asmname}/2.scaffolding/02.renaming/{asmname}.unsorted.unoriented.fa"
    output:
        "results/{asmname}/2.scaffolding/02.renaming/{asmname}.chr.unsorted.unoriented.fa"
    log:
        "results/logs/2.scaffolding/obtain_chromosomes/{asmname}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/obtain_chromosomes/{asmname}.txt"
    conda:
        "../../envs/seqkit.yaml"
    shell:
        "seqkit grep -rp \"Chr\" {input} > {output} 2> {log}"

rule rough_mashmap:
    input:
        assembly = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.chr.unsorted.unoriented.fa",
        reference = get_ref_genome,
    output:
        reference = temporary(expand("results/{{asmname}}/2.scaffolding/02.renaming/{reference}.fa", reference=config["ref_genome"])),
        out = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.chr.unsorted.unoriented.mashmap.out",
    log:
        "results/logs/2.scaffolding/rough_mashmap/{asmname}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/rough_mashmap/{asmname}.txt"
    params:
        segment = 100000,
    threads:
        len(config["ref_chr"])  #the number of chromosomes
    conda:
        "../../envs/mashmap.yaml"
    shell:
        """
        (
        ln -s $(realpath {input.reference}) {output.reference}
        samtools faidx {output.reference}
        samtools faidx {input.assembly}
        mashmap -t {threads} -s {params.segment} -r {output.reference} -q {input.assembly} -o {output.out}
        ) &> {log}
        """

rule find_wrong_orientation:
    input:
        "results/{asmname}/2.scaffolding/02.renaming/{asmname}.chr.unsorted.unoriented.mashmap.out"
    output:
        "results/{asmname}/2.scaffolding/02.renaming/{asmname}.chr.unsorted.unoriented.txt"
    log:
        "results/logs/2.scaffolding/find_wrong_orientation/{asmname}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/find_wrong_orientation/{asmname}.txt"
    shell:
        "awk 'BEGIN{{FS = OFS = \"\\t\";}} $5==\"+\"{{plus[$1]+=($4-$3); minus[$1]+=0;}} $5==\"-\"{{plus[$1]+=0; minus[$1]+=($4-$3);}} END{{for (seq in plus){{if (plus[seq]<minus[seq]){{print seq;}}}}}}' {input} > {output} 2> {log}"

rule orient_chromosomes:
    input:
        assembly = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.unsorted.unoriented.fa",
        unoriented = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.chr.unsorted.unoriented.txt",
    output:
        correct = temporary("results/{asmname}/2.scaffolding/02.renaming/{asmname}.unsorted.correct.fa.tmp"),
        unoriented = temporary("results/{asmname}/2.scaffolding/02.renaming/{asmname}.unsorted.unoriented.fa.tmp"),
        oriented = temporary("results/{asmname}/2.scaffolding/02.renaming/{asmname}.unsorted.oriented.fa.tmp"),
        combined = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.unsorted.fa",
    log:
        "results/logs/2.scaffolding/orient_chromosomes/{asmname}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/orient_chromosomes/{asmname}.txt"
    conda:
        "../../envs/seqkit.yaml"
    shell:
        """
        (
        seqkit grep -vf {input.unoriented} {input.assembly} > {output.correct}
        seqkit grep -f {input.unoriented} {input.assembly} > {output.unoriented}
        seqkit seq -rp {output.unoriented} > {output.oriented}
        cat {output.correct} {output.oriented} > {output.combined}
        ) &> {log}
        """

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
