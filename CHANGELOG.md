# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v1.8.0]
### Changed
- Updated pass/fail badges in report sample status table.
- Minor decrease to some memory directives to avoid “Process requirement exceeds available memory” errors when running in WSL.
- Default trim length set to 0 since trimming is already done by MinKNOW.
### Fixed
- The workflow no longer fails if the aliases in the sample sheet are numbers when primers are provided.
- Increase memory of `inserts` process to avoid pipeline terminated with an error exit status (137).
- Plannotate annotation infernal index error, retries and then gracefully fails for the sample, with error message provided in the wf-clone-validation-report.html.
- Fixed over-deconcatenation of assembly when there are multiple repeats. Previously the assembly could end up shorter than the provided approximate length due to other sections of sequence being removed instead of repeats.
### Added
- `min_quality` parameter to set the minimum average quality required for reads to be used in assembly (default 9). This adds more flexibility to the workflow and can improve assembly when high quality reads are available.

## [v1.7.3]
### Fixed
- Typos in "Full construct QC" section of report
- Increased mutation rate prior for bcftools call in assembly_comparison to prevent filtering of variants.
- Hidden sections weren't correctly hidden in MinKNOW
### Changed
- Linearisation table includes rows for not applicable samples.

## [v1.7.2]
This is a re-release of v1.7.1 with support for MinKNOW integration.
### Changed
- Added MinKNOW support to workflow schema.

## [v1.7.1]
### Changed
- Reconciled workflow with wf-template v5.3.4.
- Report wording to correctly reflect order of host filtering and downsampling.

## [v1.7.0]
### Added
- Support for `host_reference` and `regions_bedfile` columns in the sample sheet (to provide different references for individual samples). Sample sheet specification cannot be used in conjunction with `--host_reference` and `--regions_bedfile` parameters.
- Host reference BAM and BAI are now published in outputs.

## [v1.6.0]
### Added
- Expected reference and insert columns added to the sample status table with a tick or cross indicating if the assembly is as expected based on the provided reference.
- `--expected_identity` parameter to define the minimum expected identity percentage between provided references and the assembly (default: 99).
- `--expected_coverage` parameters for define the minimum expected coverage percentage between a provided reference and the assembly (default: 95).
- Insert QC section with coverage and BLAST Identities for Insert if an insert reference is provided.
### Changed
- Assembly reorientation attempted so that it starts at full_reference start position and matches the direction, when --full_reference is provided.
- Demo has been updated to include 12 samples.
- Emit mafs to enable dot plots to be re-used in reports generated in aliased workflows.

## [v1.5.0]
### Fixed
- Parsing of an insert that is split between the start and end of the assembly and is a reverse complement.
- Swap read count plot axis so Sample aliases are readable.
- Incorrectly running Insert QC and outputting Insert statistics when an insert was not present in the assembly.
### Changed
- Updated Medaka to v2.0.0
- The default maximum allowed mismatches in the insert primers has changed from 3 to 2.

## [v1.4.0]
### Added
- `--override_basecaller_cfg` parameter for cases where automatic basecall model detection fails or users wish to override the automatic choice.
- `--medaka_model_path` parameter to provide a custom medaka model. This is intended for users testing experimental Medaka models and will not be needed for general use. 
### Removed
- The now redundant `--basecaller_cfg` parameter as its value is now automatically detected from the input data on a per-sample basis.
### Changed
- Min and max read length determined per sample based on `approx_size`
- Emit assembly quality stats as part of final workflow outputs
- Trim length parameter can be set to 0.

## [v1.3.1]
### Added
- Parameter to control number of mismatches allowed in cutsite analysis
### Fixed
- Regression causing incorrect raw read counts shown in "Read stats" section of the report.
- INDELS now called from assembly to reference alignment (`-m1` added to bcftools mpileup)

## [v1.3.0]
### Changed
- Updated Medaka to v1.12.0.
- Updated EZCharts to v0.10.0.
### Added
- `--full_reference` parameter to accept a reference of the full construct. If provided, an additional construct QC section will be output in the report which will include reference coverage and percentage identity per sample.
- Additional per-sample output files if `--full_reference` is provided:
    - BAM: Reference aligned with the assembly in an indexed BAM.
    - Variant stats: BCF stats report with any variants found between reference and assembly.
    - Variant BCF: BCF file with all variants found between reference and assembly.
    - BAM stats: Stats report from alignment of provided reference with assembly.
- Optional linearisation efficiency section in the report, added when the `cut_site` column is supplied in the sample sheet.
- Support for `full_reference` and `insert_reference` columns in the sample sheet (to provide different references for individual samples). The MSA for insert analysis will group samples based on `insert_reference`.
- Workflow now additionally accept BAM as input by using the `--bam` parameter.
### Fixed
- Update plannotate version to v1.2.2 which fixes error that occurs when a feature contains a float. 
### Removed
- The `--medaka_model` parameter (since the appropriate Medaka model is now automatically determined from the input data). If the input data are lacking information on the model that was used to basecall them, the basecall model must be provided with `--basecaller_cfg`. Otherwise the workflow will fail.


## [v1.2.0]
### Added
- Reinstate Canu assembler alongside Flye.
- `--assembly_tool` parameter with options `canu` and `flye` (default: flye).
- `--client_fields` parameter to add extra info to the report

## [v1.1.0]
### Changed
- Ensure repetitive regions are marked on the dot plot by reducing the threshold for suppressing repeats inside exact matches.
### Added
- `--large_construct` parameter for assembly of larger constructs including Bacterial Artificial Constructs(50,000-300,000bps).

## [v1.0.0]
### Removed
- Parameters `--min_barcode` and `--max_barcode`
- Default local executor CPU and RAM limits.
### Changed
- Parameterised Flye meta option (`--non_uniform_coverage`) for non-uniform data and defaulted to false.
### Fixed
- Now handles sample aliases consisting only of numbers.

## [v0.5.3]
### Fixed
- Squashed assembly stats section downsampled plots.
### Added
- Log a warning when sample sheet approx size column is being used instead of approx size parameter.

## [v0.5.2]
### Fixed
- Deconcatenate only if the assembly is not of the approx. expected size.
- `.gbk` output file has the actual sequence for the origin
### Added
- Dotplot allowing to visualize the repetitive regions in assemblies.
- If approximate size is <=3000 set Flye min overlap to 1000.

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

