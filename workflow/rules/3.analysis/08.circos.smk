rule circos_ticks:
    output:
        "results/{asmname}/3.analysis/08.circos/{asmname}.ticks.conf",
    log:
        "results/logs/3.analysis/3.analysis/circos_ticks/{asmname}.log"
    benchmark:
        "results/benchmarks/3.analysis/circos_ticks/{asmname}.txt"
    shell:
        """
        (
        cat > {output} << EOL
show_ticks          = yes
show_tick_labels    = yes

<ticks>
skip_first_label = no
skip_last_label  = no
radius           = dims(ideogram,radius_outer)
tick_separation  = 2p
label_separation = 5p
multiplier       = 1e-6
color            = black
thickness        = 4p
size             = 20p

<tick>
spacing        = 10u
show_label     = yes
label_size     = 20p
label_offset   = 10p
format         = %d
grid           = yes
grid_color     = dgrey
grid_thickness = 1p
grid_start     = 0.5r
grid_end       = 0.999r
</tick>

</ticks>
EOL
        ) &> {log}
        """

rule circos_karyotype:
    input:
        "results/{asmname}/2.scaffolding/02.renaming/{asmname}.fa.fai"
    output:
        "results/{asmname}/3.analysis/08.circos/{asmname}.karyotype.txt",
    log:
        "results/logs/3.analysis/3.analysis/circos_karyotype/{asmname}.log"
    benchmark:
        "results/benchmarks/3.analysis/circos_karyotype/{asmname}.txt"
    shell:
        "awk 'FNR%3==0{{col=\"grey\";}} FNR%3==1{{col=\"dgrey\";}} FNR%3==2{{col=\"black\";}} {{print \"chr\", \"-\", $1, $1, \"0\", $2, col}}' {input} > {output} 2> {log}"

rule circos_configuration:
    input:
        prot_queries = expand("results/{{asmname}}/3.analysis/06.bcovblp/{query_name}.vs.{{asmname}}.items.circos", query_name=config["prot_queries"]),  ### CIRCOS ITEMS ###
        nucl_queries = expand("results/{{asmname}}/3.analysis/07.bcovbln/{query_name}.vs.{{asmname}}.fract.circos", query_name=config["nucl_queries"]),  ### CIRCOS ITEMS ###
        karyotype = "results/{asmname}/3.analysis/08.circos/{asmname}.karyotype.txt",
    output:
        "results/{asmname}/3.analysis/08.circos/{asmname}.circos.conf",
    log:
        "results/logs/3.analysis/3.analysis/circos_configuration/{asmname}.log"
    benchmark:
        "results/benchmarks/3.analysis/circos_configuration/{asmname}.txt"
    shell:
        "python3 workflow/scripts/create_circos_config.py -k {input.karyotype} -o {output} {input.prot_queries} {input.nucl_queries} &> {log}"

rule circos:
    input:
        config = "results/{asmname}/3.analysis/08.circos/{asmname}.circos.conf",
        ticks = "results/{asmname}/3.analysis/08.circos/{asmname}.ticks.conf",
    output:
        report("results/{asmname}/3.analysis/08.circos/{asmname}.circos.png", category="Circos", labels={"assembly": "{asmname}"}),
    log:
        "results/logs/3.analysis/3.analysis/circos/{asmname}.log"
    benchmark:
        "results/benchmarks/3.analysis/circos/{asmname}.txt"
    conda:
        "../../envs/circos.yaml"
    shell:
        """
        (
        cd $(dirname {output})
        circos -conf $(basename {input.config}) -outputfile $(basename {output} | rev | cut -d"." -f2- | rev)
        ) &> {log}
        """
