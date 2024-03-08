include: "01.ntjoin.smk"

rule scaffold:
    input:
        "results/{asmname}/2.scaffolding/{asmname}.vs.{ref_gen}.{minlen}.k{k}.w{w}.n1.assigned.scaffolds.fa",  #TODO: rename sequences based on mummer
    output:
        touch("results/2.scaffolding/.done")
