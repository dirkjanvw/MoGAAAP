This contact map shows the scaffolding effort done by the HapHiC tool. The
scaffolding is based on the Hi-C data and the assembly (reference-free). Please
use manual curation with e.g. Juicer to improve the scaffolding. A bash script
to create the necessary files for Juicer can be found at `results/{{ snakemake.wildcards.asmname }}/2.scaffolding/01.haphic/{{ snakemake.wildcards.asmname }}_HapHiC/04.build/juicebox.sh`.