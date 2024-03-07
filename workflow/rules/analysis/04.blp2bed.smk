rule blp2bed:
    input:
        "results/02.blast_p/{asmname}/{query_name}.vs.{asmname}.m7" 
    output: 
        "results/04.blp2bed/{asmname}/{query_name}.vs.{asmname}.bed"
    log:
        "results/logs/04.blp2bed/{asmname}/{query_name}.vs.{asmname}.log" 
    benchmark:
        "results/benchmarks/04.blp2bed/{asmname}/{query_name}.vs.{asmname}.txt" 
    conda:
        "../../envs/bedtools.yaml"
    shell:
        """
        (
        tclsh workflow/scripts/tcl_blast2bed_converter.tcl {input} {output}._temp1 1 12 13 0 BED 
        sort -k1,1 -k2,2n -k3,3n {output}._temp1 > {output}._temp2 
        bedtools merge -d 600 -i {output}._temp2 > {output} 
        perl -p -i -e 's/^ref\\|//' {output} 
        perl -p -i -e 's/\\|\t/\t/' {output} 
        ) &> {log} 
        """
