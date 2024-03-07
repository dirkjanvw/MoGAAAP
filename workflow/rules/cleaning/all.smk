include: "00.asm_bed.smk"
include: "01.blastdb.smk"
include: "02.blast_p.smk"
include: "03.blast_n.smk"
include: "04.blp2bed.smk"
include: "05.bln2bed.smk"
include: "06.bcovblp.smk"
include: "07.bcovbln.smk"

rule clean:
    input: 
        expand("results/cleaning/00.asm_bed/{asmname}/{asmname}.Asm_Len.BED", asmname = config["assemblies"]), ### BED FILE WITH SEQ LEN ###
        expand("results/cleaning/00.asm_bed/{asmname}/{asmname}.Chr_Len.BED", asmname = config["assemblies"]), ### BED FILE WITH CHR LEN ###
        expand("results/cleaning/00.asm_bed/{asmname}/{asmname}.Asm_Len.1Mb.Range", asmname = config["assemblies"]), ### ASM FILE 1MB RANGE ###
        expand("results/cleaning/00.asm_bed/{asmname}/{asmname}.Chr_Len.1Mb.Range", asmname = config["assemblies"]), ### CHR FILE 1MB RANGE ###

        expand("results/cleaning/01.blastdb/{asmname}/{asmname}.BDB", asmname = config["assemblies"]), ### ASSEMBLY BLAST DB ### 
        expand("results/cleaning/02.blast_p/{asmname}/{query_name}.vs.{asmname}.m7", asmname = config["assemblies"], query_name = config["prot_queries"]), ### BLAST PROT QUERIES ### 
        expand("results/cleaning/03.blast_n/{asmname}/{query_name}.vs.{asmname}.m7", asmname = config["assemblies"], query_name = config["nucl_queries"]), ### BLAST NUCL QUERIES ### 
        expand("results/cleaning/04.blp2bed/{asmname}/{query_name}.vs.{asmname}.bed", asmname = config["assemblies"], query_name = config["prot_queries"]), ### BED FROM PROT BLAST M7 ### 
        expand("results/cleaning/05.bln2bed/{asmname}/{query_name}.vs.{asmname}.bed", asmname = config["assemblies"], query_name = config["nucl_queries"]), ### BED FROM NUCL BLAST M7 ### 

        expand("results/cleaning/06.bcovblp/{asmname}/{query_name}.vs.{asmname}.asm.coverage", asmname = config["assemblies"], query_name = config["prot_queries"]), ### ASM COVERAGE ### 
        expand("results/cleaning/06.bcovblp/{asmname}/{query_name}.vs.{asmname}.chr.coverage", asmname = config["assemblies"], query_name = config["prot_queries"]), ### CHR COVERAGE ### 
        expand("results/cleaning/06.bcovblp/{asmname}/{query_name}.vs.{asmname}.items.circos", asmname = config["assemblies"], query_name = config["prot_queries"]), ### CIRCOS ITEMS ###

        expand("results/cleaning/07.bcovbln/{asmname}/{query_name}.vs.{asmname}.asm.coverage", asmname = config["assemblies"], query_name = config["nucl_queries"]), ### ASM COVERAGE ###
        expand("results/cleaning/07.bcovbln/{asmname}/{query_name}.vs.{asmname}.chr.coverage", asmname = config["assemblies"], query_name = config["nucl_queries"]), ### CHR COVERAGE ###
        expand("results/cleaning/07.bcovbln/{asmname}/{query_name}.vs.{asmname}.fract.circos", asmname = config["assemblies"], query_name = config["nucl_queries"]), ### CIRCOS FRACT ###
    output:
        touch("results/cleaning/.done")
