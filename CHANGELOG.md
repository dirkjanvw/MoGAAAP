All notable changes to MoGAAAP will be documented in this file.

## [0.2.1] - 2024-11-22

## Added
- Add clean GFF file: only necessary attributes and clean names (#60).

## Changed
- Remove invalid ORFs from Liftoff output (#60).

## [0.2.0] - 2024-11-07

## Added
- Add optional use of Hi-C during assembly and for visualising ntJoin contact map (#56).
- Add k-mer completeness statistics (#57).
- Use Illumina instead of HiFi for final stats if available (#57).

## [0.1.2] - 2024-10-24

## Fixed
- Also fixed OMArk, BUSCO and statistics for the assembly of a single genome (#52).
- Sometimes the mashmap environment broke due to lack of BLIS library, now added specifically to the mashmap environment (#52).
- Due to updates to numpy, BUSCO was failing, now fixed by using a newer version of BUSCO (#52).

## [0.1.1] - 2024-10-18

### Fixed
- Allow set of size 1 in config (by ignoring them effectively) (#49)

## [0.1.0] - 2024-10-14

- First release (#43)
