rule extract_organelle:
    input:
        assembly = "final_output/{asmname}.full.fa",
        asm_coverage = "results/{asmname}/2.annotation/11.bcovbln/{organelle}.vs.{asmname}.asm.coverage",
    output:
        txt = temporary("results/{asmname}/2.annotation/13.separate_genome/{asmname}.{organelle}.txt"),
        fa = "results/{asmname}/2.annotation/13.separate_genome/{asmname}.{organelle}.fa",
    log:
        "results/logs/2.annotation/extract_organelle/{asmname}/{organelle}.log",
    benchmark:
        "results/benchmarks/2.annotation/extract_organelle/{asmname}/{organelle}.txt",
    params:
        min_coverage = 0.5,  # consider hits covering more than 50% of the contig as organelle
    conda:
        "../../envs/seqkit.yaml"
    shell:
        """
        (
        awk -v mincov={params.min_coverage} '{{if ($7 >= mincov) print $1}}' {input.asm_coverage} > {output.txt}
        seqkit grep -f {output.txt} {input.assembly} > {output.fa}
        ) &> {log}
        """

rule extract_nuclear:
    input:
        assembly = "final_output/{asmname}.full.fa",
        organelles = expand("results/{{asmname}}/2.annotation/13.separate_genome/{{asmname}}.{organelle}.fa", organelle=config["organellar"]),
    output:
        txt = temporary("results/{asmname}/2.annotation/13.separate_genome/{asmname}.nuclear.txt"),
        fa = "results/{asmname}/2.annotation/13.separate_genome/{asmname}.nuclear.fa",
    log:
        "results/logs/2.annotation/extract_nuclear/{asmname}.log",
    benchmark:
        "results/benchmarks/2.annotation/extract_nuclear/{asmname}.txt",
    conda:
        "../../envs/seqkit.yaml"
    shell:
        """
        (
        awk '/^>/{{print substr($1,2)}}' {input.organelles} > {output.txt}
        seqkit grep -v -f {output.txt} {input.assembly} > {output.fa}
        ) &> {log}
        """
