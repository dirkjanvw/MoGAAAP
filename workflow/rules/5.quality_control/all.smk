rule qc:
    input:
        ""
    output:
        touch("results/quality_control.done")
