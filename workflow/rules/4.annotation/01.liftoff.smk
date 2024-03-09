rule liftoff:
    input:
        ref_annotation = config["ref_annotation"],
        ref_genome = get_ref_genome,
        assembly = "results/{asmname}/2.scaffolding/03.renaming/{asmname}.fa",
    output:
        refgff = temporary("results/{asmname}/4.annotation/01.liftoff/reference.gff"),
        gff = protected("results/{asmname}/4.annotation/01.liftoff/liftoff.gff"),
        polished = protected("results/{asmname}/4.annotation/01.liftoff/liftoff.gff_polished"),
    log:
        "results/logs/4.annotation/liftoff/{asmname}.log"
    benchmark:
        "results/benchmarks/4.annotation/liftoff/{asmname}.txt"
    threads:
        10
    conda:
        "../../envs/liftoff.yaml"
    shell:
        """
        (
        ln -s $(realpath {input.gff}) {output.refgff}
        liftoff -p {threads} -copies -cds -polish -u $(dirname {output.polished})/unmapped_features.txt -dir $(dirname {output.polished})/intermediate_files -o {output.gff} -g {output.refgff} {input.query} {input.ref}
        ) &> {log}
        """
