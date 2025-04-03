rule bcovblp:
    input:
        asm_len_bed  = "results/{asmname}/2.annotation/04.asm_bed/{asmname}.Asm_Len.BED",
        chr_1Mb_blk  = "results/{asmname}/2.annotation/04.asm_bed/{asmname}.Chr_Len.1Mb.Range",
        bed_file     = "results/{asmname}/2.annotation/08.blp2bed/{query_name}.vs.{asmname}.bed"
    output:
        asm_coverage = "results/{asmname}/2.annotation/10.bcovblp/{query_name}.vs.{asmname}.asm.coverage",
        chr_coverage = "results/{asmname}/2.annotation/10.bcovblp/{query_name}.vs.{asmname}.chr.coverage",
        circos_file  = "results/{asmname}/2.annotation/10.bcovblp/{query_name}.vs.{asmname}.items.circos"
    log:
        "results/logs/2.annotation/bcovblp/{asmname}/{query_name}.vs.{asmname}.log"
    benchmark:
        "results/benchmarks/2.annotation/bcovblp/{asmname}/{query_name}.vs.{asmname}.txt"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        """
        (
        bedtools coverage -a {input.asm_len_bed} -b {input.bed_file} > {output.asm_coverage}
        bedtools coverage -a {input.chr_1Mb_blk} -b {input.bed_file} > {output.chr_coverage}
        cut -f 1,2,3,4 {output.chr_coverage} > {output.circos_file}
        ) &> {log}
        """
