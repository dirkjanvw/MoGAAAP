#!/bin/bash
# This is a wrapper script around the MoGAAAP Snakemake pipeline
# It is intended as a convenience script to run the pipeline
#
# Usage: mogaaap <command> [options]
#
# Commands:
#   configure   Configure the configuration files using command line arguments
#   validate    Validate the configuration files and existence of SIF files
#   run         Run the MoGAAAP pipeline


# Print usage
usage() {
    cat <<EOF

=== MoGAAAP wrapper script ===

Usage: $(basename "$0") <command> [options]

Commands:
  configure   Configure the pipeline
                NB: This only works when all samples have the same reference genome
    Options:
      -s, --samples        Sample sheet TSV file
                             (required)
      -f, --reference      Reference genome FASTA file
                             (required)
      -g, --annotation     Reference annotation GFF3 file
                             (required)
      -m, --helixer_model  Helixer model file
                             (required)
      -x, --gxdb           GXDB database location
                             (required)
      -b, --odb            ODB name for BUSCO
                             (required)
      -k, --kraken_db      Kraken2 database location
                             (required)
      -a, --OMA-db         OMA database location
                             (required)
      -t, --telomere_motif Telomere motif
                             default: CCCTAAA
      -o, --output         Output configuration YAML file
                             default: config.yaml

  validate    Validate the configuration and pipeline
    Options:
      -v, --verbose  Print verbose output
                       default: false

  run         Run the MoGAAAP pipeline
    Options:
      -t, --target   Target to run (all, assemble, scaffold, analyse, annotate,
                           qa)
                       default: all
      -r, --report   Generate a report
                       default: false
      -n, --dryrun  Perform a dry-run
                       default: false
      -c, --cores    Maximum allowed number of cores
                       default: 10
      -m, --memory   Maximum allowed memory (in GB)
                       default: 10

Global options:
  -h, --help  Print this help message and exit

EOF
    exit 1
}


# Handle configure subcommand
handle_configure() {
    # Parse options
    telomere_motif=CCCTAAA
    output=config.yaml
    while getopts ":s:f:g:m:x:b:k:a:t:o:" opt; do
        case $opt in
            s) samples="$OPTARG" ;;
            f) reference="$OPTARG" ;;
            g) annotation="$OPTARG" ;;
            m) helixer_model="$OPTARG" ;;
            x) gxdb="$OPTARG" ;;
            b) odb="$OPTARG" ;;
            k) kraken="$OPTARG" ;;
            a) OMAdb="$OPTARG" ;;
            t) telomere_motif="$OPTARG" ;;
            o) output="$OPTARG" ;;
            \?) echo "Error: Invalid option -$OPTARG" >&2; usage ;;
            :) echo "Error: Option -$OPTARG requires an argument" >&2; usage ;;
        esac
    done

    # Check if required options are given
    if [[ -z "${samples}" || -z "${reference}" || -z "${annotation}" || -z "${helixer_model}" || -z "${gxdb}" || -z "${odb}" || -z "${kraken}" || -z "${OMAdb}" ]]; then
        echo "Error: Missing required options"
        usage
    fi

    # Check if files exist or can be created
    fileerror=false
    if [[ ! -f "${samples}" ]]; then
        echo "Error: Sample sheet file not found: ${samples}"
        fileerror=true
    fi
    if [[ ! -f "${reference}" ]]; then
        echo "Error: Reference genome file not found: ${reference}"
        fileerror=true
    fi
    if [[ ! -f "${annotation}" ]]; then
        echo "Error: Reference annotation file not found: ${annotation}"
        fileerror=true
    fi
    if [[ ! -f "${helixer_model}" ]]; then
        echo "Error: Helixer model file not found: ${helixer_model}"
        fileerror=true
    fi
    if [[ ! -d "${gxdb}" ]]; then
        echo "Error: GXDB directory not found: ${gxdb}"
        fileerror=true
    fi
    if [[ ! -d "${kraken}" ]]; then
        echo "Error: Kraken2 database directory not found: ${kraken}"
        fileerror=true
    fi
    if [[ ! -d "${OMAdb}" ]]; then
        echo "Error: OMA database directory not found: ${OMAdb}"
        fileerror=true
    fi
    if [[ -f "${output}" ]]; then
        echo "Error: Output file already exists: ${output}"
        fileerror=true
    fi
    if ${fileerror}; then
        exit 1
    fi

    # Get reference name from first sample
    refname=$(awk 'BEGIN{FS = OFS = "\t";} FNR==1{for (i=1;i<=NF;i++) if ($i == "referenceId") {col=i; break;}} FNR==2{print $col; exit;}' ${samples})

    # Configure the pipeline
    echo "Writing configuration to ${output}"
    echo "samples: $(realpath ${samples})" > ${output}
    echo "assembler: hifiasm" >> ${output}
    echo "scaffolder: ntjoin" >> ${output}
    echo "min_contig_len: 10000" >> ${output}
    echo "ntjoin_k: 52" >> ${output}
    echo "ntjoin_w: 16000" >> ${output}
    echo "reference_genomes:" >> ${output}
    echo "    ${refname}:" >> ${output}
    echo "        genome: $(realpath ${reference})" >> ${output}
    echo "        annotation: $(realpath ${annotation})" >> ${output}
    echo "        chromosomes:" >> ${output}
    awk '/^>/{printf "            %d: \"%s\"\n", ++c, substr($0,2);}' ${reference} >> ${output}
    echo "#prot_queries: #optional" >> ${output}
    echo "#nucl_queries: #optional" >> ${output}
    echo "#organellar: #optional" >> ${output}
    echo "telomere_motif: \"${telomere_motif}\"" >> ${output}
    echo "helixer_model: $(realpath ${helixer_model})" >> ${output}
    echo "helixer_max_gene_length: 64152" >> ${output}
    echo "k_qa: 21" >> ${output}
    echo "gxdb: $(realpath ${gxdb})" >> ${output}
    echo "pantools_grouping: 3" >> ${output}
    echo "odb: ${odb}" >> ${output}
    echo "kraken2_nt: $(realpath ${kraken})" >> ${output}
    echo "OMAdb: $(realpath ${OMAdb})" >> ${output}
    echo "set:" >> ${output}
    echo "    all:" >> ${output}
    for id in $(cut -f 1 ${samples} | tail -n +2); do
        echo "        - ${id}" >> ${output}
    done
    echo "jvm: \"-Xmx100g -Xms100g\"" >> ${output}
    echo "tmpdir: /dev/shm" >> ${output}
    echo "nucmer_maxgap: 1000" >> ${output}
    echo "nucmer_minmatch: 1000" >> ${output}
}


# Handle validate subcommand
handle_validate() {
    # Parse options
    verbose=false
    while getopts ":v" opt; do
        case $opt in
            v) verbose=true ;;
            \?) echo "Error: Invalid option -$OPTARG" >&2; usage ;;
        esac
    done

    # Validate the pipeline
    echo "Validating pipeline"
    valid=true
    for def in $(ls workflow/singularity/*/*def); do
        if [[ -f "${def%def}sif" ]]; then
            if $verbose; then
                echo "Found SIF file for $(basename $(dirname ${def})): ${def%def}sif"
            fi
        else
            echo "Error: SIF file not found for $(basename $(dirname ${def})): ${def%def}sif"
            valid=false
        fi
    done

    # Exit if pipeline is not valid
    if ! $valid; then
        exit 1
    fi

    # Validate the configuration by running Snakemake with --dryrun
    echo "Validating configuration"
    echo snakemake ${target} \
      --printshellcmds \
      --cores ${cores} \
      --dryrun \
      --use-conda --use-singularity \
      --resources gbmem=${memory} helixer=1 pantools=1
}


# Handle run subcommand
handle_run() {
    # Parse options
    target=all
    dryrun=
    report=false
    cores=10
    memory=10
    while getopts ":nrt:c:m:" opt; do
        case $opt in
            n) dryrun=true ;;
            r) report=true ;;
            t) target="$OPTARG" ;;
            c) cores="$OPTARG" ;;
            m) memory="$OPTARG" ;;
            \?) echo "Error: Invalid option -$OPTARG" >&2; usage ;;
            :) echo "Error: Option -$OPTARG requires an argument" >&2; usage ;;
        esac
    done

    # Check if target is valid
    if [[ ! "$target" =~ ^(all|assemble|scaffold|analyse|annotate|qa)$ ]]; then
        echo "Error: Invalid target: $target"
        usage
    fi

    # Run the pipeline
    echo snakemake ${target} \
      --printshellcmds \
      --cores ${cores} \
      ${dryrun:+--dryrun} \
      --use-conda --use-singularity \
      --resources gbmem=${memory} helixer=1 pantools=1

    # Create report if requested
    if $report; then
        echo snakemake ${target} --report report.html
    fi
}


# Check if any arguments are given
if [[ $# -eq 0 ]]; then
    usage
fi


# Handle global options
while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            ;;
        *)
            break
            ;;
    esac
done


# Handle subcommands
subcommand="$1"
shift
case "$subcommand" in
    configure)
        handle_configure "$@"
        ;;
    validate)
        handle_validate "$@"
        ;;
    run)
        handle_run "$@"
        ;;
    *)
        echo "Error: Unknown command: $subcommand"
        usage
        ;;
esac
