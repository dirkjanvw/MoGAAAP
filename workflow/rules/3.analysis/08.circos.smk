rule circos_ticks:
    output:
        "results/{asmname}/3.analysis/08.circos/ticks.conf",
    log:
        "results/logs/3.analysis/circos_ticks/{asmname}.log"
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
        "results/logs/3.analysis/circos_karyotype/{asmname}.log"
    benchmark:
        "results/benchmarks/3.analysis/circos_karyotype/{asmname}.txt"
    shell:
        "awk '$1!~/Chr[0-9]+/{{next;}} c%3==0{{col=\"black\";}} c%3==1{{col=\"dgrey\";}} c%3==2{{col=\"grey\";}} {{print \"chr\", \"-\", $1, $1, \"0\", $2, col; c++;}}' {input} > {output} 2> {log}"

rule circos_configuration:
    input:
        counts = expand("results/{{asmname}}/3.analysis/06.bcovblp/{query_name}.vs.{{asmname}}.items.circos", query_name=config["prot_queries"]),  ### CIRCOS ITEMS ###
        fractions = expand("results/{{asmname}}/3.analysis/07.bcovbln/{query_name}.vs.{{asmname}}.fract.circos", query_name=config["nucl_queries"]),  ### CIRCOS ITEMS ###
        karyotype = "results/{asmname}/3.analysis/08.circos/{asmname}.karyotype.txt",
        ticks = "results/{asmname}/3.analysis/08.circos/ticks.conf",
    output:
        conf = "results/{asmname}/3.analysis/08.circos/{asmname}.circos.conf",
        overview = report("results/{asmname}/3.analysis/08.circos/{asmname}.circos.tsv", category="Circos", labels={"file": "overview", "assembly": "{asmname}"}),
    log:
        "results/logs/3.analysis/circos_configuration/{asmname}.log"
    benchmark:
        "results/benchmarks/3.analysis/circos_configuration/{asmname}.txt"
    shell:
        """
        (
        printf "filename\\tfiletype\\tmin\\tmax\\n" > {output.overview}
        for file in {input.counts}; do
            printf "${{file}}\\tcount\\t0\\t100\\n" >> {output.overview}
        done
        for file in {input.fractions}; do
            printf "${{file}}\\tfraction\\t0\\t1.0\\n" >> {output.overview}
        done
        printf "{input.karyotype}\\tkaryotype\\tNA\\tNA\\n" >> {output.overview}
        ln -s $(realpath {input.counts} {input.fractions}) $(dirname {output.conf})/
        SCRIPT=$(realpath workflow/scripts/create_circos_config.py)
        cd $(dirname {output.conf})
        python3 $SCRIPT -k $(basename {input.karyotype}) -o $(basename {output.conf}) $(echo {input.counts} {input.fractions} | awk '{{for (i=1;i<=NF;i++){{n=split($i,a,"/"); print a[n];}}}}')
        ) &> {log}
        """

rule circos:
    input:
        "results/{asmname}/3.analysis/08.circos/{asmname}.circos.conf",
    output:
        png = report("results/{asmname}/3.analysis/08.circos/{asmname}.circos.png", category="Circos", labels={"file": "plot", "assembly": "{asmname}"}),
        svg = "results/{asmname}/3.analysis/08.circos/{asmname}.circos.svg",
    log:
        "results/logs/3.analysis/circos/{asmname}.log"
    benchmark:
        "results/benchmarks/3.analysis/circos/{asmname}.txt"
    conda:
        "../../envs/circos.yaml"
    shell:
        """
        (
        cd $(dirname {output.png})
        circos -conf $(basename {input}) -outputfile $(basename {output.png} | rev | cut -d"." -f2- | rev)
        ) &> {log}
        """
