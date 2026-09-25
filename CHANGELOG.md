# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## UNRELEASED
### Added
- [GCGI-1666](https://jira.oicr.on.ca/browse/GCGI-1666) - CRAM support: the reference from the selected `genomeVersion` is passed to every task, including `-ref_genome` for AMBER and COBALT so CRAM inputs can be decoded
- Optional `input_amber_directory` and `input_cobalt_directory`; when set the corresponding task is skipped and PURPLE reads from the pre-computed directory
- New `genomeVersion` option `grch38_hmf` added for the HMF GRCh38_masked_exclusions_alts_hlas reference
### Changed
- Renamed required workflow inputs: `tumour_bam`/`tumour_bai`/`normal_bam`/`normal_bai` are now `tumour`/`tumour_index`/`normal`/`normal_index`, since they accept bam or cram
- `runPURPLE` inputs `amber_directory`/`cobalt_directory` split into `amber_zip`/`cobalt_zip` (in-workflow task output) and `amber_dir`/`cobalt_dir` (pre-computed directory)
- `amber_directory` and `cobalt_directory` outputs are now optional, as those tasks may be skipped
- `genomeVersion` now defaults to `grch38_hmf`
- `filterSMALL` defaults now point at the grch38-hmf reference and modules, in line with the `genomeVersion` default
- Updated regression tests, README and commands.txt for the renamed and added parameters

## [1.4.0] - 2026-02-26
### Changed
- [GBS-6821](https://jira.oicr.on.ca/browse/GBS-6821) - continue to switch to new versions of tools
- Upgrade purple, cobalt and amber to match oncoanalyser 2.3.0 pipeline 
### Fixed
- fixed parameters which were incompatible with new versions

## [1.3.1] - 2026-02-03
### Changed
- [GRD-1021](https://jira.oicr.on.ca/browse/GRD-1021) - continue to switch to new versions of tools
- disabled LINX
- doSV flag is still in the workflow but filtering is conditioned on the existance of SV file input
### Fixed
- fixed parameters which were incompatible with new versions

## [1.3.0] - 2026-01-10
### Changed
- switched to hmftools/1.2 with tools of the same versions as used by Heartwig pipeline at the time of comparison to GSI pipeline

## [1.2.4] - 2026-01-10
### Added
- Added ability to configure Java Heap for tasks

## [1.2.3] - 2025-12-18
### Added
- Added more reference options, noAlt and ncbi references 

## [1.2.2] - 2025-09-32
### Fixed
- Fixed a bug, incorrect filtering of the small variants input (snv vcf from mutect2)
- [GRD-993](https://jira.oicr.on.ca/browse/GRD-993)
### Changed
- Modified wdl metadata and README in accordance with the correction
- Modified wdl metadata/dependencies to indicate additional tools used in the workflow

## [1.2.1] - 2025-05-26
- Re-deployment to enable labels for optional outputs
- [GRD-948](https://jira.oicr.on.ca/browse/GRD-948)

## [1.2.0] - 2024-06-25
### Added
- [GRD-797](https://jira.oicr.on.ca/browse/GRD-797) - Add vidarr labels to outputs (changes to medata only)

## [1.1.3] - 2024-02-16
### Added
- Added alternate solution

## [1.1.2] - 2023-11-06
### Changed
- Parameter changes, regression test and documents

## [1.1.1] - 2023-11-03
### Added
- Min diploid tumor ratio count

## [1.0] - 2023-09-18
### Changed
- Aggregated updates

## [0.0] - 2023-02-01
### Added
- A brand-new workflow.
