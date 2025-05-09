rule bln2bed:
    input:
        "results/{asmname}/2.annotation/07.blast_n/{query_name}.vs.{asmname}.m7"
    output:
        "results/{asmname}/2.annotation/09.bln2bed/{query_name}.vs.{asmname}.bed"
    log:
        "results/logs/2.annotation/bln2bed/{asmname}/{query_name}.vs.{asmname}.log"
    benchmark:
        "results/benchmarks/2.annotation/bln2bed/{asmname}/{query_name}.vs.{asmname}.txt"
    conda:
        "../../envs/bedtools.yaml"
    shell:
        """
        (
        tclsh workflow/scripts/tcl_blast2bed_converter.tcl {input} {output}._temp1 1 12 13 0 BED
        sort -k1,1 -k2,2n -k3,3n {output}._temp1 > {output}._temp2
        bedtools merge -d 100 -i {output}._temp2 > {output}
        ### REMOVE (PERL /FIND/REPLACE/) "ref|" and "|" SYMBOLS IN OUTPUT FILE ###
        perl -p -i -e 's/^ref\\|//' {output}
        perl -p -i -e 's/\\|\t/\t/' {output}
        ) &> {log}
        """
