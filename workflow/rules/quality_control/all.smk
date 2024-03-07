rule qc:
    input:
        ""
    output:
        touch("results/5.quality_control/.done")
