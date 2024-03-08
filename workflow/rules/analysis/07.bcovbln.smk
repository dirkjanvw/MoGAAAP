rule bcovbln:
    input: 
        asm_len_bed = "results/{asmname}/3.analysis/00.asm_bed/{asmname}.Asm_Len.BED",
        chr_1Mb_blk = "results/{asmname}/3.analysis/00.asm_bed/{asmname}.Chr_Len.1Mb.Range",
        bed_file = "results/{asmname}/3.analysis/05.bln2bed/{query_name}.vs.{asmname}.bed",
    output:
        asm_coverage = "results/{asmname}/3.analysis/07.bcovbln/{query_name}.vs.{asmname}.asm.coverage",
        chr_coverage = "results/{asmname}/3.analysis/07.bcovbln/{query_name}.vs.{asmname}.chr.coverage",
        circos_file = "results/{asmname}/3.analysis/07.bcovbln/{query_name}.vs.{asmname}.fract.circos",
    log:
        "results/logs/bcovbln/{asmname}/{query_name}.vs.{asmname}.log"
    benchmark:
        "results/benchmarks/bcovbln/{asmname}/{query_name}.vs.{asmname}.txt"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        """
        ( 
        bedtools coverage -a {input.asm_len_bed} -b {input.bed_file} > {output.asm_coverage} 
        bedtools coverage -a {input.chr_1Mb_blk} -b {input.bed_file} > {output.chr_coverage} 
        cut -f 1,2,3,7 {output.chr_coverage} > {output.circos_file} 
        ) &> {log} 
        """
