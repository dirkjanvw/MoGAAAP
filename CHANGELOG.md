All notable changes to MoGAAAP will be documented in this file.

## [Unreleased]

### Added
- Add option to use the pipeline for already existing assemblies (#83).

### Changed
- Update sample sheet to allow using the pipeline for already existing assemblies (#83).

## [0.2.6 - 2025-02-21]

### Added
- Add `test_data` directory with full end-to-end user test instructions (#86).

## [0.2.5 - 2025-02-10]

### Added
- Add option to use the pipeline without custom built singularity containers (#79).
- Add option to use the pipeline for already existing assemblies (#83).

### Changed
- No longer require sudo rights for setup (#81).
- Update sample sheet to allow using the pipeline for already existing assemblies (#83).
- Added some more information to README about the pipeline (#82).

## [0.2.4 - 2025-01-20]

### Added
- Messages on start, success and error (#73).
- RagTag is now available for scaffolding as alternative to ntJoin (#75).

## [0.2.3 - 2024-12-18]

### Added
- Creation of extra large MUMmerplot (#70).

## [0.2.2 - 2024-12-11]

### Fixed
- Fix PanTools upset plot creation for accessions with 2 haplotypes (#65).
- Fix the BUSCO plot when there are 2 haplotypes (#66).
- Prevent duplicate chromosome names (#67).

### Added
- More clearly mark in report which WGS data was used for which statistics (#64).

## [0.2.1] - 2024-11-22

### Added
- Add clean GFF file: only necessary attributes and clean names (#60).

### Changed
- Remove invalid ORFs from Liftoff output (#60).

## [0.2.0] - 2024-11-07

### Added
- Add optional use of Hi-C during assembly and for visualising ntJoin contact map (#56).
- Add k-mer completeness statistics (#57).
- Use Illumina instead of HiFi for final stats if available (#57).

## [0.1.2] - 2024-10-24

### Fixed
- Also fixed OMArk, BUSCO and statistics for the assembly of a single genome (#52).
- Sometimes the mashmap environment broke due to lack of BLIS library, now added specifically to the mashmap environment (#52).
- Due to updates to numpy, BUSCO was failing, now fixed by using a newer version of BUSCO (#52).

## [0.1.1] - 2024-10-18

### Fixed
- Allow set of size 1 in config (by ignoring them effectively) (#49)

## [0.1.0] - 2024-10-14

- First release (#43)
