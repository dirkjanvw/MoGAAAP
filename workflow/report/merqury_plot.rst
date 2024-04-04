This plot is a so-called spectra-CN plot. What is shows is the distribution of
k-mers in a given WGS dataset. These k-mers are then coloured according to their
occurrence in the assembly. The x-axis (*kmer_multiplicity*) is a proxy for the
coverage of the WGS dataset for the organism. The y-axis (*Count*) is the number
of k-mers for a given multiplicity.

Since the colour indicates the occurrence in the assembly, the plot can be used
to assess the completeness of the assembly. Read-only k-mers should ideally only
be present at very low multiplicity and not at the heterozygous and homozygous
peaks. The heterozygous peak should be made up of k-mers that are present 1x in
the assembly and the homozygous peak should be made up of k-mers that are
present 2x in the assembly.

It is important to realise the type of WGS data used for the plot:

- **HiFi**: These high-quality long-reads are used for assembling the genome and
  should therefore show a "perfect" picture described above.
- **ONT**: As these long-reads are generally of lower quality than HiFi reads,
  ONT reads are generally only informative of a good quality genome when they
  have a very high coverage. This allows for the visual separation of the
  sequencing artefact, heterozygous and homozygous peaks.
- **Illumina**: These reads are not used in this pipeline for the assembly at
  any point and are therefore a good independent measure of the quality of the
  assembly. A peak of read-only k-mers at the heterozygous peak indicates a poor
  assembly or a mismatch between the assembly and the WGS data (mislabelled
  sample or contamination). This should be accompanied by a high number of
  k-mers at multiplicity 0 that are present (>0) in the assembly.