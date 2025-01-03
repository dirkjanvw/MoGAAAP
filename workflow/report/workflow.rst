This report has been created by the `MoGAAAP
<https://github.com/dirkjanvw/MoGAAAP>`_ and summarises the results
of all assembled, scaffolded and annotated genomes. Please bear in mind that
beyond minimal contig length, no filtering has been applied to the assemblies
or quick-and-dirty annotations. The goal of this report is to provide an
overview of the quality of the assemblies and annotations, so that the user may
decide on further steps.

This report contains the following pages:

- **Workflow**: The current page. On the right, an overview of all steps in the
  pipeline is shown. Each step is clickable and will take you to a detailed
  technical overview of that step.
- **Statistics**: The pipeline statistics: when was what step run and how long
  it took. This page is not always accurate nor complete.
- **About**: Technical details of this report; only relevant to Snakemake
  developers.
- **Result**: All possible categories of the results:

  - *Analysis*: Results for user-defined sequence searches.
  - *Circos*: A circos plot per assembly, to be used as template.
  - *Collinearity*: Collinearity between sets of assemblies as generated by
    `ntSynt`.
  - *Contamination*: Several contamination assessments for each assembly.
  - *Gene completeness*: Gene completeness assessments by BUSCO and OMArk for
    each assembly.
  - *General statistics*: General statistics for sets of assemblies.
  - *K-mer completeness*: K-mer completeness assessments by BUSCO and OMArk for
    each assembly.
  - *Mapping*: Mapping of Illumina reads to the assembly (if available).
  - *MUMmerplot*: A MUMmerplot per assembly, aligned to the defined reference
    genome.
  - *Pangrowth*: K-mer based pangenome analysis of sets of assemblies (if more
    than two assemblies are present).
  - *PanTools*: Gene based pangenome analysis of sets of assemblies (if more
    than two assemblies are present).
  - *Phylogeny*: Several phylogenetic trees for sets of assemblies.
