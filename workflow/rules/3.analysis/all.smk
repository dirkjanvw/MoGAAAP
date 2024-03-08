include: "00.asm_bed.smk"
include: "01.blastdb.smk"
include: "02.blast_p.smk"
include: "03.blast_n.smk"
include: "04.blp2bed.smk"
include: "05.bln2bed.smk"
include: "06.bcovblp.smk"
include: "07.bcovbln.smk"

rule analyse:
    input:
        expand("results/{asmname}/3.analysis/00.asm_bed/{asmname}.Asm_Len.BED", asmname = get_all_accessions()), ### BED FILE WITH SEQ LEN ###
        expand("results/{asmname}/3.analysis/00.asm_bed/{asmname}.Chr_Len.BED", asmname = get_all_accessions()), ### BED FILE WITH CHR LEN ###
        expand("results/{asmname}/3.analysis/00.asm_bed/{asmname}.Asm_Len.1Mb.Range", asmname = get_all_accessions()), ### ASM FILE 1MB RANGE ###
        expand("results/{asmname}/3.analysis/00.asm_bed/{asmname}.Chr_Len.1Mb.Range", asmname = get_all_accessions()), ### CHR FILE 1MB RANGE ###

        expand("results/{asmname}/3.analysis/02.blast_p/{query_name}.vs.{asmname}.m7", asmname = get_all_accessions(), query_name = config["prot_queries"]), ### BLAST PROT QUERIES ###
        expand("results/{asmname}/3.analysis/03.blast_n/{query_name}.vs.{asmname}.m7", asmname = get_all_accessions(), query_name = config["nucl_queries"]), ### BLAST NUCL QUERIES ###
        expand("results/{asmname}/3.analysis/04.blp2bed/{query_name}.vs.{asmname}.bed", asmname = get_all_accessions(), query_name = config["prot_queries"]), ### BED FROM PROT BLAST M7 ###
        expand("results/{asmname}/3.analysis/05.bln2bed/{query_name}.vs.{asmname}.bed", asmname = get_all_accessions(), query_name = config["nucl_queries"]), ### BED FROM NUCL BLAST M7 ###

        expand("results/{asmname}/3.analysis/06.bcovblp/{query_name}.vs.{asmname}.asm.coverage", asmname = get_all_accessions(), query_name = config["prot_queries"]), ### ASM COVERAGE ###
        expand("results/{asmname}/3.analysis/06.bcovblp/{query_name}.vs.{asmname}.chr.coverage", asmname = get_all_accessions(), query_name = config["prot_queries"]), ### CHR COVERAGE ###
        expand("results/{asmname}/3.analysis/06.bcovblp/{query_name}.vs.{asmname}.items.circos", asmname = get_all_accessions(), query_name = config["prot_queries"]), ### CIRCOS ITEMS ###

        expand("results/{asmname}/3.analysis/07.bcovbln/{query_name}.vs.{asmname}.asm.coverage", asmname = get_all_accessions(), query_name = config["nucl_queries"]), ### ASM COVERAGE ###
        expand("results/{asmname}/3.analysis/07.bcovbln/{query_name}.vs.{asmname}.chr.coverage", asmname = get_all_accessions(), query_name = config["nucl_queries"]), ### CHR COVERAGE ###
        expand("results/{asmname}/3.analysis/07.bcovbln/{query_name}.vs.{asmname}.fract.circos", asmname = get_all_accessions(), query_name = config["nucl_queries"]), ### CIRCOS FRACT ###

        #TODO: add circos plot of the above data
    output:
        touch("results/3.analysis/.done")
