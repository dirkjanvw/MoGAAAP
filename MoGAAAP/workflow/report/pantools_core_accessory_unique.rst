This growth plot shows the pangenome growth curve of the pangenome defined by
"{{ snakemake.wildcards.asmset }}". Homology is based on the protein sequences
of the transcripts encoded by the genes. The curves are drawn by random sampling
a given number of genomes (x-axis) and then counting how many homology groups
are core, accessory or unique. The accessory and unique genome together make up
the dispensable genome. (Unique means unique to a genome and accessory means
neither core nor unique.)

The core genome is expected to stabilise at some point, depending on how many
genomes are sampled. When this core curve flattens off, these homology groups
make up the core genome of the pangenome. The dispensable genome is expected to
grow as more genomes are sampled. In case of a open pangenome, this growth is
expected to continue indefinitely; a closed pangenome will eventually reach a
plateau.