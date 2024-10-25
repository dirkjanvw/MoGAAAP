This contact map shows the scaffolding effort done by the HapHiC tool. The
scaffolding is based on the Hi-C data and the assembly (reference-free). In the
contact map, the axes represent the groups of scaffolds that were scaffolded
together. The color intensity represents the number of Hi-C contacts between
the groups of scaffolds. These contacts are used to determine the physical
distance between sequences in the assembly.

Please use manual curation with e.g. Juicer to improve the scaffolding. A bash
script to create the necessary files for Juicer can be found at `results/{{ snakemake.wildcards.asmname }}/2.scaffolding/01.haphic/{{ snakemake.wildcards.asmname }}_HapHiC/04.build/juicebox.sh`
Ideally, a diagonal should be seen from the left bottom to the top right, with
the centromeres, telomeres and other repeats as off-diagonal signals. In case
other off-diagonal signals are present, these indicate mis-joins during the
scaffolding process (or misassemblies in the assembly).