rule circos:
    input:
        "results/{asmname}/3.analysis/00.asm_bed/{asmname}.Asm_Len.BED",  ### BED FILE WITH SEQ LEN ###
        "results/{asmname}/3.analysis/00.asm_bed/{asmname}.Chr_Len.BED",  ### BED FILE WITH CHR LEN ###
        "results/{asmname}/3.analysis/00.asm_bed/{asmname}.Asm_Len.1Mb.Range",  ### ASM FILE 1MB RANGE ###
        "results/{asmname}/3.analysis/00.asm_bed/{asmname}.Chr_Len.1Mb.Range",  ### CHR FILE 1MB RANGE ###

        expand("results/{{asmname}}/3.analysis/06.bcovblp/{query_name}.vs.{{asmname}}.asm.coverage", query_name=config["prot_queries"]),  ### ASM COVERAGE ###
        expand("results/{{asmname}}/3.analysis/06.bcovblp/{query_name}.vs.{{asmname}}.chr.coverage", query_name=config["prot_queries"]),  ### CHR COVERAGE ###
        expand("results/{{asmname}}/3.analysis/06.bcovblp/{query_name}.vs.{{asmname}}.items.circos", query_name=config["prot_queries"]),  ### CIRCOS ITEMS ###

        expand("results/{{asmname}}/3.analysis/07.bcovbln/{query_name}.vs.{{asmname}}.asm.coverage", query_name=config["nucl_queries"]),  ### ASM COVERAGE ###
        expand("results/{{asmname}}/3.analysis/07.bcovbln/{query_name}.vs.{{asmname}}.chr.coverage", query_name=config["nucl_queries"]),  ### CHR COVERAGE ###
    output:
        "results/{asmname}/3.analysis/08.circos/{asmname}.circos.png",
    log:
        "results/logs/3.analysis/circos/{asmname}.log"
    benchmark:
        "results/benchmarks/3.analysis/circos/{asmname}.txt"
    conda:
        "../../envs/circos.yaml"
    shell:
        "circos ??? &> {log}"