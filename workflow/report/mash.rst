This plot shows a mash-based distance between all assemblies of a given set. As
this plot is created by `mash`, it only considers a subset of all k-mers in the
assembly. It is therefore good to compare this plot to
{{ snakemake.wildcards.asmname }}.csv.mash.pdf_ which is created from all k-mers in
the assembly.

Ideally, the relationship between the assemblies should correspond to the known
phylogenetic relationship between the organisms. If this is not the case, it may
point to mislabeling of samples or contamination.