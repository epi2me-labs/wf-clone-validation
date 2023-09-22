# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unrelease]
### Fixed
- `.gbk` output file has the actual sequence for the origin

## [v0.5.1]
### Added
- Additionally output plannotate annotations as a GenBank (`.gbk`) file.
### Removed
- Remove plasmid length column from plannotate feature table.
### Changed
- Strand column of plannotate feature to use `+` and `-` notation.
- Default `--primers` parameter is now set to null.
### Fixed
- If `--host_filter` provided the fastcat stats will be included in the report.

## [v0.5.0]
### Changed
- The report has been updated and re-ordered to improve usability.
- Default basecaller cfg is now `dna_r10.4.1_e8.2_400bps_sup@v4.2.0`.
- Docker will use an ARM platform image on appropriate devices.
### Fixed
- Updated basecaller cfg model options.
- Workflow will still output report if there are no assemblies.

## [v0.4.0]
### Added
- Documentation updated to include workflow steps.
- Full plasmid assembly mean quality table in report.
- Output a fastq of the final assembly.
- Insert reference, if provided, will now be used to variant call insert consensus with bcftools.

### Removed
- Unused packages from the container.

### Changed
- Enum choices are enumerated in the `--help` output
- Enum choices are enumerated as part of the error message when a user has selected an invalid choice
- Bumped minimum required Nextflow version to 22.10.8
- Updated GitHub issue templates to force capture of more information.
- Reference parameter changed to `--insert_reference`.
- Updated example command displayed when running `--help`
- Parameter`--approx_size_sheet` no longer accepted, instead use sample sheet with optional additional column `approx_size`. 
- Any sample aliases that contain spaces will be replaced with underscores.

### Fixed
- Replaced `--threads` option in fastqingress with hardcoded values to remove warning about undefined `param.threads`
- Annotation output bed file has correct notation for strand.

## [v0.3.1]
### Added
- Configuration for running demo data in AWS

## [v0.3.0]
### Added
- Flye replaces canu as the assembler tool.

## [v0.2.13]
### Changed
- Updated to Oxford Nanopore Technologies PLC. Public License.

### Fixed
- Amended raw QC stats to show data before filtering by assembly_size parameter.

## [v0.2.12]
### Fixed
- Bug where the workflow wouldn't run properly when `--approx_size_sheet` was used.

### Changed
- Now uses new `fastq_ingress` implementation.

## [v0.2.11]
### Fixed
- Provide medaka model for each assembly to fix bug.

## [v0.2.10]
### Fixed
- Replace spaces with tabs in medaka model TSV to fix bug.

## [v0.2.9]
### Fixed
- Medaka models added to container

## [v0.2.8]
### Changed
- `--basecall_cfg` is now used to determine suitable Medaka model, alternatively provide the name of a model with `--medaka_model` to override automatic selection.

## [v0.2.7]
### Changed
- Updated description in manifest

## [v0.2.6]
### Removed
- `-profile conda` is no longer supported, users should use `-profile standard` (Docker) or `-profile singularity` instead

### Added
- `nextflow run epi2me-labs/wf-clone-validation --version` will now print the workflow version number and exit

## [v0.2.5]
### Fixed
- Filter host step not outputting approx_size.

### Updated
- Use groovy script to ping after workflow has run.

## [v0.2.4]
### Added
- Error handling for no annotations found for an assembly.
- Windows parameter so Canu can run on windows

### Fixed
- Plannotate dictionary keys can contain any characters.
- Sanitize fastq intermittent null object error.

## [v0.2.3]
### Changed
- Change params.threads to task.cpus

## [v0.2.2]
### Changed
- Fastqingress metadata map.
- Sample status now collected from tuples.

## [v0.2.1]
### Changed
- Plannotate read/write database requirement fix
- approx_size_sheet param instead of sample_sheet
- Set out_dir option type to ensure output is written to correct directory on Windows

## [v0.2.0]
### Changed
- Better help text on CLI.
- Fix issue with S3 file inputs.

### Updated
- Plannotate to version v1.2.0

## [v0.1.9]
### Added
- Param for fast option in Canu assembly

## [v0.1.8]
### Changed
- New docs format

## [v0.1.7]
### Fixed
- Sample sheet encoding
- Min max barcodes integer types

### Changed
- Moved bioinformatics from report to seperate processes

### Added
- Ability to define approx_size of sequence per sample in sample_sheet
- Insert length to table
- Output annotation bed files per sample

## [v0.1.6]
### Changed
- Update schema for epi2melabs compatibility

### Fixed
- Make use of the canu_useGrid parameter

## [v0.1.5]
### Added
- Singularity profile to config.
- Ping telemetry file.
- Handle more fastq input directory structures.

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
* First release

