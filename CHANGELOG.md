# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [v0.1.5]
### Added
- Singularity profile to config.
### Fixed
- db_directory description and explained in README
- db_directory param updated to match s3 folder name
### Changed
- Use downsampled samples for polish assembly step.

## [v0.1.4]
### Added
- Option to add suffix to HTML report name.
- Error message if fastq input file evaluates to null.

## [v0.1.3]
### Fixed
- Default Primers parameter txt to tsv.

## [v0.1.2]
### Added
- Fastcat stats plots in tabs for pass and failed samples.
- Version and parameter tables.
- Per barcode number of reads.
- Insert sequences output.
- MSA of inserted sequences.
### Fixed
- Order samples lexicographically.
### Changed
- Use Canu for assembly instead of Flye.
- Trim input sequences.

## [v0.1.1]
### Fixed
- Corrected number of input channels for host_reference process.
- Remove duplicate output files.
- Help message parameters reflect config.

## [v0.1.0]
### Added
- Plannotate for plasmid annotation and visualization.
- Per sample pass or fail error message in CSV.
- Plasmid annotation feature table output CSV.
### Changed
- Updated project to use latest practices from wf-template.
### Fixed
- Incorrect specification of conda environment file in Nextflow config.

## [v0.0.4]
### Added
- Fix report naming to be consistent with other projects

## [v0.0.3]
### Added
- Optional --prefix flag for naming outputs

## [v0.0.2]
### Added
- --no-reconcile flag for a simpler and quicker overall pipeline

### Changed
- simplified assembly outputs, to only emit the final polished assembly


## [v0.0.1]

First release
