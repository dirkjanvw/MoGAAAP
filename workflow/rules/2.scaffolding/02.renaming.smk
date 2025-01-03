rule create_renaming_table:
    input:
        agp = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.agp",
        mxdot = "results/{asmname}/2.scaffolding/01.ntjoin/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.mx.dot",
    output:
        "results/{asmname}/2.scaffolding/02.renaming/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.conversion.tsv"
    log:
        "results/logs/2.scaffolding/create_renaming_table/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.log"
    benchmark:
        "results/benchmarks/2.scaffolding/create_renaming_table/{asmname}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.txt"
    params:
        num_chr = lambda wildcards: len(config["reference_genomes"][get_reference_id(wildcards.asmname)]["chromosomes"]),
    shell: # I found the following awk commands to get a quick list of what chromosomes belong together (based on https://github.com/bcgsc/ntJoin/issues/63):
        """
        (
        awk 'BEGIN{{OFS = "'\\''";}} {{print $1,$6;}}' {input.agp} | \
            awk -vasm="{wildcards.asmname}" -vref="{wildcards.reference}" 'BEGIN{{FS = "'\\''"; OFS = "\\t";}} FNR==NR{{chr[$2]=$1; next;}} $1~"^"ref{{printf "%s\\t",$2;}} $1~"^"asm{{printf "%s\\n",chr[$2];}}' - {input.mxdot} | \
            sort | \
            uniq -c | \
            sort -nr | \
            awk -vnum_chr="{params.num_chr}" 'BEGIN{{OFS = "\\t";}} FNR<=num_chr && !a[$2] {{print $2,$3; a[$2]=1;}}' > {output}
        ) &> {log}
        """

rule visualise_ntjoin_renaming:
    input:
        table = lambda wildcards: expand("results/{{asmname}}/2.scaffolding/02.renaming/{{asmname}}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.conversion.tsv",
            reference=get_reference_id(wildcards.asmname),
            minlen=config["min_contig_len"],
            k=config["ntjoin_k"],
            w=config["ntjoin_w"],
        ),
    output:
        report("results/{asmname}/2.scaffolding/02.renaming/{asmname}.html",
            category="Hi-C",
            labels={"assembly": "{asmname}",
                    "stage": "scaffolds",
                    "algorithm": "ntJoin (conversion table)"}
        ),
    log:
        "results/logs/2.scaffolding/visualise_ntjoin_renaming/{asmname}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/visualise_ntjoin_renaming/{asmname}.txt"
    conda:
        "../../envs/csvtotable.yaml"
    shell:
        """
        (
        awk 'BEGIN{{FS = OFS = \"\\t\"; print \"chromosome name\",\"scaffold name\";}} {{print;}}' {input.table} > {input.table}.tmp
        csvtotable -d $'\\t' -o {input.table}.tmp {output}
        rm {input.table}.tmp
        ) &> {log}
        """

rule renaming_scaffolds:
    input:
        all = lambda wildcards: expand("results/{{asmname}}/2.scaffolding/01.ntjoin/{{asmname}}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.fa",
            reference=get_reference_id(wildcards.asmname),
            minlen=config["min_contig_len"],
            k=config["ntjoin_k"],
            w=config["ntjoin_w"],
        ),
        table = lambda wildcards: expand("results/{{asmname}}/2.scaffolding/02.renaming/{{asmname}}.vs.{reference}.min{minlen}.k{k}.w{w}.n2.all.scaffolds.conversion.tsv",
            reference=get_reference_id(wildcards.asmname),
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
        chroms = lambda wildcards: config["reference_genomes"][get_reference_id(wildcards.asmname)]["chromosomes"],
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
        reference = lambda wildcards: config["reference_genomes"][get_reference_id(wildcards.asmname)]["genome"],
    output:
        reference = temporary("results/{asmname}/2.scaffolding/02.renaming/reference.fa"),
        out = "results/{asmname}/2.scaffolding/02.renaming/{asmname}.chr.unsorted.unoriented.mashmap.out",
    log:
        "results/logs/2.scaffolding/rough_mashmap/{asmname}.log"
    benchmark:
        "results/benchmarks/2.scaffolding/rough_mashmap/{asmname}.txt"
    params:
        segment = 100000,
    threads:
        lambda wildcards: len(config["reference_genomes"][get_reference_id(wildcards.asmname)]["chromosomes"]),  #the number of chromosomes
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
        "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa",
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
