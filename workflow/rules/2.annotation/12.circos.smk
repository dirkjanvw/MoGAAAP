rule circos_ticks:
    output:
        "results/{asmname}/2.annotation/12.circos/ticks.conf",
    log:
        "results/logs/2.annotation/circos_ticks/{asmname}.log"
    benchmark:
        "results/benchmarks/2.annotation/circos_ticks/{asmname}.txt"
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
        "final_output/{asmname}.full.fa.fai"
    output:
        "results/{asmname}/2.annotation/12.circos/{asmname}.karyotype.txt",
    log:
        "results/logs/2.annotation/circos_karyotype/{asmname}.log"
    benchmark:
        "results/benchmarks/2.annotation/circos_karyotype/{asmname}.txt"
    shell: #circos has a maximum of 200 sequences it can show, so limit to either "Chr[0-9]+" containing sequences or the first 100 sequences
        "awk '$1!~/Chr[0-9]+/ && FNR > 100{{next;}} c%3==0{{col=\"black\";}} c%3==1{{col=\"dgrey\";}} c%3==2{{col=\"grey\";}} {{print \"chr\", \"-\", $1, $1, \"0\", $2, col; c++;}}' {input} > {output} 2> {log}"

def get_circos_files(wildcards):
    prot_files = []
    if "prot_queries" in config:
        for query in config["prot_queries"]:
            prot_files.append(f"results/{wildcards.asmname}/2.annotation/10.bcovblp/{query}.vs.{wildcards.asmname}.items.circos")

    nucl_files = []
    if "nucl_queries" in config:
        for query in config["nucl_queries"]:
            nucl_files.append(f"results/{wildcards.asmname}/2.annotation/11.bcovbln/{query}.vs.{wildcards.asmname}.fract.circos")

    # If you want to add organellar files, uncomment the following lines
    #organellar_files = []
    #if "organellar" in config:
    #    for query in config["organellar"]:
    #        organellar_files.append(f"results/{wildcards.asmname}/2.annotation/11.bcovbln/{query}.vs.{wildcards.asmname}.fract.circos")

    return prot_files + nucl_files #+ organellar_files

rule circos_configuration:
    input:
        files = get_circos_files,  ### CIRCOS ITEMS ###
        karyotype = "results/{asmname}/2.annotation/12.circos/{asmname}.karyotype.txt",
        ticks = "results/{asmname}/2.annotation/12.circos/ticks.conf",
    output:
        conf = "results/{asmname}/2.annotation/12.circos/{asmname}.circos.conf",
        overview = "results/{asmname}/2.annotation/12.circos/{asmname}.circos.tsv",
    log:
        "results/logs/2.annotation/circos_configuration/{asmname}.log"
    benchmark:
        "results/benchmarks/2.annotation/circos_configuration/{asmname}.txt"
    shell:
        """
        (
        printf "filename\\tfiletype\\tmin\\tmax\\n" > {output.overview}
        for file in {input.files}; do
            if [[ "${{file}}" =~ \\.fract\\.circos$ ]]; then
                printf "${{file}}\\tfraction\\t0\\t1.0\\n" >> {output.overview}
            else
                printf "${{file}}\\tcount\\t0\\t100\\n" >> {output.overview}
            fi
        done
        printf "{input.karyotype}\\tkaryotype\\tNA\\tNA\\n" >> {output.overview}
        cp {input.files} $(dirname {output.conf})/
        SCRIPT=$(realpath workflow/scripts/create_circos_config.py)
        cd $(dirname {output.conf})
        python3 $SCRIPT -k $(basename {input.karyotype}) -o $(basename {output.conf}) $(echo {input.files} | awk '{{for (i=1;i<=NF;i++){{n=split($i,a,"/"); print a[n];}}}}')
        ) &> {log}
        """

rule visualise_circos_configuration:
    input:
        "results/{asmname}/2.annotation/12.circos/{asmname}.circos.tsv"
    output:
        report("results/{asmname}/2.annotation/12.circos/{asmname}.circos.html",
            category="Circos",
            caption="../../report/circos_overview.rst",
            labels={"file": "overview", "assembly": "{asmname}"}),
    log:
        "results/logs/2.annotation/visualise_circos_configuration/{asmname}.log"
    benchmark:
        "results/benchmarks/2.annotation/visualise_circos_configuration/{asmname}.txt"
    conda:
        "../../envs/csvtotable.yaml"
    shell:
        "csvtotable -d $'\\t' {input} {output} &> {log}"

rule circos:
    input:
        "results/{asmname}/2.annotation/12.circos/{asmname}.circos.conf",
    output:
        png = report("results/{asmname}/2.annotation/12.circos/{asmname}.circos.png",
            category="Circos",
            caption="../../report/circos_plot.rst",
            labels={"file": "plot", "assembly": "{asmname}"}),
        svg = "results/{asmname}/2.annotation/12.circos/{asmname}.circos.svg",
    log:
        "results/logs/2.annotation/circos/{asmname}.log"
    benchmark:
        "results/benchmarks/2.annotation/circos/{asmname}.txt"
    conda:
        "../../envs/circos.yaml"
    shell:
        """
        (
        cd $(dirname {output.png})
        circos -conf $(basename {input}) -outputfile $(basename {output.png} | rev | cut -d"." -f2- | rev)
        ) &> {log}
        """
