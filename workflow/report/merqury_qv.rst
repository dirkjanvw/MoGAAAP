This file contains the QV scores for each sequence in the scaffolded assembly as
calculated by `merqury`. The value is a
:math:`QV = -10 \log_{10}(\text{error rate})` score, where an error means: a
k-mer in the assembly which cannot be found in the WGS data.

It is important to realise the type of WGS data used to calculate the QV scores:

- **HiFi**: These high-quality long reads were used to assembly the genome and
  will therefore cause a relatively high QV score as few errors are expected.
- **ONT**: Even though these long reads were used to assembly the genome if
  provided, their quality is lower than the HiFi reads and will therefore cause
  a lower QV score.
- **Illumina**: These reads are not used in this pipeline for the assembly at
  any point and are therefore a good independent measure of the quality of the
  assembly. A low QV score here indicates a poor assembly or a mismatch between
  the assembly and the WGS data (mislabelled sample or contamination).

For the overall QV score, please see the {{ snakemake.wildcards.asmname }}.pdf_
file.