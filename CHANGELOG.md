# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.3.1] - 2026-02-03
### Changed
- [GRD-1021](https://jira.oicr.on.ca/browse/GRD-1021) - continue to switch to new versions of tools
- disabled LINX
- fixed parameters which were incompatible with new versions
- doSV flag is still in the workflow but filtering is conditioned on the existance of SV file input

## [1.3.0] - 2026-01-16
### Changed
- [GRD-1021](https://jira.oicr.on.ca/browse/GRD-1021)
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
