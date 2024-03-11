rule blp2bed:
    input:
        "results/{asmname}/3.analysis/02.blast_p/{query_name}.vs.{asmname}.m7"
    output:
        "results/{asmname}/3.analysis/04.blp2bed/{query_name}.vs.{asmname}.bed"
    log:
        "results/logs/3.analysis/blp2bed/{asmname}/{query_name}.vs.{asmname}.log"
    benchmark:
        "results/benchmarks/3.analysis/blp2bed/{asmname}/{query_name}.vs.{asmname}.txt"
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
