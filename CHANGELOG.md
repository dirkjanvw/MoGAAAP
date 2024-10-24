All notable changes to MoGAAAP will be documented in this file.

## [UNRELEASED]

## Fixed
- Also fixed OMArk, BUSCO and statistics for the assembly of a single genome (#52).
- Sometimes the mashmap environment broke due to lack of BLIS library, now added specifically to the mashmap environment (#52).
- Due to updates to numpy, BUSCO was failing, now fixed by using a newer version of BUSCO (#52).

## [0.1.1] - 2024-10-18

### Fixed
- Allow set of size 1 in config (by ignoring them effectively) (#49)

## [0.1.0] - 2024-10-14

- First release (#43)
