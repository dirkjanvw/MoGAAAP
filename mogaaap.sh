#!/bin/bash
# This is a wrapper script around the MoGAAAP Snakemake pipeline
# It is intended as a convenience script to run the pipeline
#
# Usage: mogaaap [options]


# Print usage
usage() {
    cat <<EOF

=== MoGAAAP wrapper script ===

Usage: $(basename "$0") [options]

Options:
  -h, --help           Print this help message and exit

  [General options]
  -y, --config         Configuration YAML file (required)

  [Configuration options (writes to --config); all required]
  --generate-config    Generate a configuration YAML file and exit
  -s, --samples        Sample sheet TSV file
  --reference-fasta    Reference genome FASTA file
  --reference-gff3     Reference annotation GFF3 file
  --mitochondrion      Mitochondrial genome FASTA file
  --chloroplast        Chloroplast genome FASTA file (optional)
  --helixer-model      Helixer model file
  --gxdb               GXDB database location
  --odb                ODB name for BUSCO
  --kraken-db          Kraken2 database location
  --OMA-db             OMA database location
  --telomere-motif     Telomere motif (default: CCCTAAA)
  --use-custom-singularity
                       Use custom singularity containers (default: true)
                       NB: Setting this to false will reduce the output

  [Pipeline options]
  --run                Run the pipeline
  -t, --target         Target to run (all, assemble, scaffold, analyse,
                            annotate, qa) (default: all)
  -r, --report         Generate a report
  -n, --dryrun         Perform a dry-run
  -c, --cores          Maximum allowed number of cores (default: 10)
  -m, --memory         Maximum allowed memory (in GB) (default: 500)

EOF
    exit 1
}


# Parse options
config=

generate_config=false
samples=
reference=
annotation=
mitochondrion=
chloroplast=
helixer_model=
gxdb=
odb=
kraken=
OMAdb=
telomere_motif=CCCTAAA
use_custom_singularity=

run=false
target=all
report=false
dryrun=false
cores=10
memory=500

while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            ;;
        -s|--samples)
            samples="$2"
            shift 2
            ;;
        -y|--config)
            config="$2"
            shift 2
            ;;
        --generate-config)
            generate_config=true
            shift
            ;;
        --reference-fasta)
            reference="$2"
            shift 2
            ;;
        --reference-gff3)
            annotation="$2"
            shift 2
            ;;
        --mitochondrion)
            mitochondrion="$2"
            shift 2
            ;;
        --chloroplast)
            chloroplast="$2"
            shift 2
            ;;
        --helixer-model)
            helixer_model="$2"
            shift 2
            ;;
        --gxdb)
            gxdb="$2"
            shift 2
            ;;
        --odb)
            odb="$2"
            shift 2
            ;;
        --kraken-db)
            kraken="$2"
            shift 2
            ;;
        --OMA-db)
            OMAdb="$2"
            shift 2
            ;;
        --telomere-motif)
            telomere_motif="$2"
            shift 2
            ;;
        --run)
            run=true
            shift
            ;;
        -t|--target)
            target="$2"
            shift 2
            ;;
        -r|--report)
            report=true
            shift
            ;;
        -n|--dryrun)
            dryrun=true
            shift
            ;;
        -c|--cores)
            cores="$2"
            shift 2
            ;;
        -m|--memory)
            memory="$2"
            shift 2
            ;;
        --use-custom-singularity)
            use_custom_singularity="$2"
            shift 2
            ;;
        *)
            echo "Error: Unknown option: $1"
            usage
            ;;
    esac
done


# Validate options
incomplete=false
if [[ -z "${config}" ]]; then
    echo "Error: Missing required option: --config"
    incomplete=true
fi
if ${generate_config} && ${run}; then
    echo "Error: --generate-config cannot be combined with --run"
    incomplete=true
fi
if ! ${generate_config} && ! ${run}; then
    echo "Error: Either --generate-config or --run must be specified"
    incomplete=true
fi
if ${generate_config}; then
    if [[ -z "${samples}" ]]; then
        echo "Error: Missing required option for --generate-config: --samples"
        incomplete=true
    fi
    if [[ -z "${reference}" ]]; then
        echo "Error: Missing required option for --generate-config: --reference-fasta"
        incomplete=true
    fi
    if [[ -z "${annotation}" ]]; then
        echo "Error: Missing required option for --generate-config: --reference-gff3"
        incomplete=true
    fi
    if [[ -z "${mitochondrion}" ]]; then
        echo "Error: Missing required option for --generate-config: --mitochondrion"
        incomplete=true
    fi
    if [[ -z "${helixer_model}" ]]; then
        echo "Error: Missing required option for --generate-config: --helixer-model"
        incomplete=true
    fi
    if [[ -z "${gxdb}" ]]; then
        echo "Error: Missing required option for --generate-config: --gxdb"
        incomplete=true
    fi
    if [[ -z "${odb}" ]]; then
        echo "Error: Missing required option for --generate-config: --odb"
        incomplete=true
    fi
    if [[ -z "${kraken}" ]]; then
        echo "Error: Missing required option for --generate-config: --kraken-db"
        incomplete=true
    fi
    if [[ -z "${OMAdb}" ]]; then
        echo "Error: Missing required option for --generate-config: --OMA-db"
        incomplete=true
    fi
    if [[ -z "${use_custom_singularity}" ]]; then
        use_custom_singularity=true
    fi
    if [[ ! "${use_custom_singularity}" =~ ^(true|false)$ ]]; then
        echo "Error: Invalid option for --use-custom-singularity: ${use_custom_singularity} (only true/false allowed)"
        incomplete=true
    fi
    if ${report}; then
        echo "Error: --report cannot be combined with --generate-config"
        incomplete=true
    fi
    if ${dryrun}; then
        echo "Error: --dryrun cannot be combined with --generate-config"
        incomplete=true
    fi
fi
if ${run}; then
    if [[ ! -f "${config}" ]]; then
        echo "Error: Configuration file not found: ${config}"
        incomplete=true
    fi
    if [[ -n "${samples}" ]]; then
        echo "Error: --samples cannot be combined with --run"
        incomplete=true
    fi
    if [[ -n "${reference}" ]]; then
        echo "Error: --reference-fasta cannot be combined with --run"
        incomplete=true
    fi
    if [[ -n "${annotation}" ]]; then
        echo "Error: --reference-gff3 cannot be combined with --run"
        incomplete=true
    fi
    if [[ -n "${mitochondrion}" ]]; then
        echo "Error: --mitochondrion cannot be combined with --run"
        incomplete=true
    fi
    if [[ -n "${chloroplast}" ]]; then
        echo "Error: --chloroplast cannot be combined with --run"
        incomplete=true
    fi
    if [[ -n "${helixer_model}" ]]; then
        echo "Error: --helixer-model cannot be combined with --run"
        incomplete=true
    fi
    if [[ -n "${gxdb}" ]]; then
        echo "Error: --gxdb cannot be combined with --run"
        incomplete=true
    fi
    if [[ -n "${odb}" ]]; then
        echo "Error: --odb cannot be combined with --run"
        incomplete=true
    fi
    if [[ -n "${kraken}" ]]; then
        echo "Error: --kraken-db cannot be combined with --run"
        incomplete=true
    fi
    if [[ -n "${OMAdb}" ]]; then
        echo "Error: --OMA-db cannot be combined with --run"
        incomplete=true
    fi
    if [[ -n "${use_custom_singularity}" ]]; then
        echo "Error: --use-custom-singularity cannot be combined with --run"
        incomplete=true
    fi
    if [[ ! "${target}" =~ ^(all|assemble|scaffold|analyse|annotate|qa)$ ]]; then
        echo "Error: Invalid target: ${target}"
        incomplete=true
    fi
fi
if ${incomplete}; then
    exit 1
fi


# Generate configuration file if requested
if ${generate_config}; then
    echo "Generating configuration file: ${config}"

    # Check existence files
    critical_error=false
    if [[ ! -f "${samples}" ]]; then
        echo "Error: Sample sheet not found: ${samples}"
        critical_error=true
    fi
    if [[ ! -f "${reference}" ]]; then
        echo "Error: Reference genome not found: ${reference}"
        critical_error=true
    fi
    if [[ ! -f "${annotation}" ]]; then
        echo "Error: Reference annotation not found: ${annotation}"
        critical_error=true
    fi
    if [[ ! -f "${mitochondrion}" ]]; then
        echo "Error: Mitochondrial genome not found: ${mitochondrion}"
        critical_error=true
    fi
    if [[ -n "${chloroplast}" ]] && [[ ! -f "${chloroplast}" ]]; then
        echo "Error: Chloroplast genome not found: ${chloroplast}"
        critical_error=true
    fi
    if [[ ! -f "${helixer_model}" ]]; then
        echo "Error: Helixer model not found: ${helixer_model}"
        critical_error=true
    fi
    if [[ ! -d "${gxdb}" ]]; then
        echo "Error: GXDB database not found: ${gxdb}"
        critical_error=true
    fi
    if [[ ! -d "${kraken}" ]]; then
        echo "Error: Kraken2 database not found: ${kraken}"
        critical_error=true
    fi
    if [[ ! -f "${OMAdb}" ]]; then
        echo "Error: OMA database not found: ${OMAdb}"
        critical_error=true
    fi
    if ${critical_error}; then
        exit 1
    fi

    # Get reference name from first sample
    refname=$(awk 'BEGIN{FS = OFS = "\t";} FNR==1{for (i=1;i<=NF;i++) if ($i == "referenceId") {col=i; break;}} FNR==2{print $col; exit;}' ${samples})

    # Create telomere motif file for blastn
    telomere_file=$(mktemp)
    echo ${telomere_motif} | awk 'BEGIN{print ">telomere";} {for (i=1;i<=100;i++){printf "%s",$1;} printf "\n";}' > ${telomere_file}

    cat <<EOF > "${config}"
# MoGAAAP configuration file (auto-generated)
# Please edit this file to customize the pipeline

# This file contains the full configuration for the pipeline in YAML format.
# Make sure to fill out at least everything that is not labelled as optional.
# All strings that start with "resources" should be relative or absolute paths
# on your machine that are readable by Snakemake.

# All samples should be defined in a TSV file. In this file, a unique
# accessionId, path to HiFi reads, (optional) path to ONT reads, (optional)
# paths to Illumina reads, species name, NCBI taxonomy ID, and reference genome
# ID should be defined.
samples: $(realpath ${samples})

# Whether custom singularity containers can be used. If set to false, not all
# output will be available.
custom_singularity: ${use_custom_singularity}

# Assembly settings
assembler: "hifiasm" #supported assemblers: hifiasm, verkko

# Scaffolding settings
scaffolder: "ntjoin" #supported scaffolders: ntjoin, ragtag

# All filtering and scaffolding parameters
min_contig_len: 10000
ntjoin_k: 52 #required for ntjoin
ntjoin_w: 16000 #required for ntjoin

# All reference genomes and chromosomes. The reference genome ID in the samples
# sheet should match the key in this section.
reference_genomes:
    ${refname}:
        genome: $(realpath ${reference})
        annotation: $(realpath ${annotation})
        chromosomes:
EOF
    # Add chromosomes to configuration
    awk '/^>/{printf "            %d: \"%s\"\n", ++c, substr($1,2);}' "${reference}" >> "${config}"
    # And continue
    cat <<EOF >> "${config}"

# All protein and nucleotide queries that need to be run for the assembly
# analysis module. The key should be the name of the query and the value should
# be the path to the query. At least one protein or nucleotide query is
# required for the module to run (organellar sequences count as nucleotide
# queries).
# Contigs with a more than 50% of its length covered by an organellar sequence
# will be separated from the scaffolded assembly output (not for QC, though).
# For the protein queries, the query should be a single exon domain as it
# is searched in the genome and not in the proteome
#prot_queries: #optional
nucl_queries: #optional but recommended
    telomere: $(realpath ${telomere_file}) #for blastn: should be a 100x repeat of the telomere motif to identify telomeres anywhere in the genome
organellar:
    mitochondrion: $(realpath ${mitochondrion})
EOF
    # Add chloroplast if available
    if [[ -n "${chloroplast}" ]]; then
        echo "    chloroplast: $(realpath ${chloroplast})" >> "${config}"
    fi
    # And continue
    cat <<EOF >> "${config}"
telomere_motif: "${telomere_motif}" #for seqtk: it only identifies the motif at the end of a chromosome

# Annotation settings
helixer_model: $(realpath ${helixer_model})
helixer_max_gene_length: 64152

# Quality assessment settings
k_qa: 21 #optimal k-mer size for quality assessment
gxdb: $(realpath ${gxdb}) #path to gxdb for fcs-gx
pantools_grouping: 3 #grouping relaxation setting for pantools
odb: ${odb} #ODB database for busco
kraken2_nt: $(realpath ${kraken}) #kraken2 database
OMAdb: $(realpath ${OMAdb}) #OMA database for orthologs

# A section to define groups of accessions for plotting together in report. This
# list is order sensitive and will only be used for QC.
set:
    all:
EOF
    # Add samples to configuration
    for id in $(cut -f 1 ${samples} | tail -n +2); do
        echo "        - ${id}" >> "${config}"
    done
    # And continue
    cat <<EOF >> "${config}"

# General settings; only change if you know what you are doing
jvm: "-Xmx100g -Xms100g" #java virtual machine settings
tmpdir: "/dev/shm/" #temporary directory for fast IO
nucmer_maxgap: 1000 #maximum gap size for nucmer
nucmer_minmatch: 1000 #minimum match size for nucmer
EOF
    exit 0

    # Log success
    echo "Configuration file generated: ${config}"
    echo "Please edit this file to customize the pipeline"
fi


# Run the pipeline if requested
if ${run}; then
    problems=false

    # Check existence of Snakemake
    if ! command -v snakemake &> /dev/null; then
        echo "Error: Snakemake not found"
        problems=true
    fi

    # Log start
    if ${dryrun}; then
        echo "Performing a dry-run"
    else
        echo "Running the pipeline"
        dryrun=
    fi

    # Check existence of Singularity or Apptainer
    if ! command -v singularity &> /dev/null; then
        if ! command -v apptainer &> /dev/null; then
            echo "Error: Singularity nor Apptainer found"
            problems=true
        fi
    fi

    # Check existence of configuration file
    if [[ ! -f "${config}" ]]; then
        echo "Error: Configuration file not found: ${config}"
        problems=true
    fi

    # If custom_singularity in config file has been set to true, check for presence of SIF files in workflow/singularity directory
    if grep -q "custom_singularity: true" ${config}; then
        for def in workflow/singularity/*/*.def; do
            if [[ ! -f "${def%def}sif" ]]; then
                echo "Error: Custom Singularity container not found: ${def%def}sif"
                problems=true
            fi
        done
    fi

    # Exit if any problems were found
    if ${problems}; then
        exit 1
    fi

    # Run the pipeline
    snakemake ${target} \
      --configfile ${config} \
      --printshellcmds \
      --cores ${cores} \
      ${dryrun:+--dryrun} \
      --use-conda --use-singularity \
      --resources gbmem=${memory}

    # Create report if requested
    if ${report}; then
        snakemake ${target} \
          ${dryrun:+--dryrun} \
          --report report.html
    fi
fi
