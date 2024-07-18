rule asm_bed:
    input:
        "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa",
    output:
        asm_bed = "results/{asmname}/3.analysis/00.asm_bed/{asmname}.Asm_Len.BED",
        chr_bed = "results/{asmname}/3.analysis/00.asm_bed/{asmname}.Chr_Len.BED",
        asm_1Mb = "results/{asmname}/3.analysis/00.asm_bed/{asmname}.Asm_Len.1Mb.Range",
        chr_1Mb = "results/{asmname}/3.analysis/00.asm_bed/{asmname}.Chr_Len.1Mb.Range",
    log:
        "results/logs/3.analysis/asm_bed/{asmname}.log"
    benchmark:
        "results/benchmarks/3.analysis/asm_bed/{asmname}.txt"
    params:
        num_chr = lambda wildcards: len(config["reference_genomes"][get_reference_id(wildcards.asmname)]["chromosomes"])
    conda:
        "../../envs/bioawk.yaml"
    shell:
        """
        (
        bioawk -c fastx '{{ print $name, 0, length($seq) }}' {input} > {output.asm_bed}
        head -n {params.num_chr} {output.asm_bed} > {output.chr_bed}
        python2 workflow/scripts/make_subrange_from_bed.py 1000000 < {output.asm_bed} > {output.asm_1Mb}
        python2 workflow/scripts/make_subrange_from_bed.py 1000000 < {output.chr_bed} > {output.chr_1Mb}
        ) &> {log}
        """
