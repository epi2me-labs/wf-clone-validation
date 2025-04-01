### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| fastq | string | FASTQ files to use in the analysis. | This accepts one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| bam | string | BAM or unaligned BAM (uBAM) files to use in the analysis. | This accepts one of three cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files; (iii) the path to a directory containing one level of sub-directories which in turn contain BAM files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| approx_size | integer | Approximate size of the plasmid in base pairs. This can also be defined per sample: see sample_sheet parameter. |  | 7000 |
| assm_coverage | integer | Fold coverage for use per assembly | This is the coverage that will be used to subsample reads to use for the assembly. | 60 |
| primers | string | TSV File containing primers used to find inserts. If left empty then inserts will not be searched for. | Specify one or more primer sets which will be used to find the sequence inserted in the construct. This file should be in .tsv format containing columns [primer_name,  5' primer, 3' primer] with no header. An example `primers.tsv` for pRham/T7 is available in the demo data for the workflow. |  |
| analyse_unclassified | boolean | Analyse unclassified reads from input directory. By default the workflow will not process reads in the unclassified directory. | If selected and if the input is a multiplex directory the workflow will also process the unclassified directory. | False |


### Reference Genome Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| insert_reference | string | Optional file containing insert reference sequence which will be used for comparison with consensus insert in the report. | Providing a reference sequence can be useful as a QC on the base-level resolution of the the reconstructed insert sequences. Users can specify different insert references for individual samples using the sample sheet and including an `insert_reference` column. This cannot be used in conjunction with `--insert_reference`. |  |
| full_reference | string | Optional FASTA file containing the reference sequence of the full plasmid. This will be used for comparison with the assembled construct. | Providing a reference sequence can be useful as a quality check on the base-level resolution of the reconstructed sequence, the reference is not used to generate the assembly. Users can specify different full references for individual samples using the sample sheet and including a `full_reference` column. This cannot be used in conjunction with `--full_reference`. |  |
| host_reference | string | A host reference genome FASTA file. Read which map to this reference are discarded and not used for the assembly.  Users can specify different host references for individual samples using the sample sheet and including a `host_reference` column. This cannot be used in conjunction with `--host_reference`. |  |  |
| regions_bedfile | string | If a host_reference supplied, add an optional BED file to provide host reference regions that will be masked during filtering.  Users can specify different BED files for individual samples using the sample sheet and including a `regions_bedfile` column. This cannot be used in conjunction with `--regions_bedfile`. |  |  |


### Sample Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_sheet | string | A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a directory containing sub-directories with FASTQ files. An optional column `approx_size` can be added to provide size estimates for each sample. When not provided, the `--approx_size` parameter will be used for all samples. | The sample sheet is a CSV file with, minimally, columns named `barcode` and `alias`. Extra columns are allowed. A `type` column is required for certain workflows and should have the following values; `test_sample`, `positive_control`, `negative_control`, `no_template_control`. The workflow can use a number of optional columns: `approx_size` provides optional size estimates for each sample, `cut_site` can be added to provide a cut site as a sequence which will be used to provide a linearisation efficiency section in the report, `full_reference` and `insert_reference` allow the use of per-sample references when providing full/relative paths (with respect to the workflow launch directory) to the respective reference files. |  |
| sample | string | A single sample name for non-multiplexed data. Permissible if passing a single .fastq(.gz) file or directory of .fastq(.gz) files. |  |  |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all user-facing files. |  | output |
| prefix | string | The prefix attached to each of the output filenames. |  |  |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| override_basecaller_cfg | string | Override auto-detected basecaller model that processed the signal data; used to select an appropriate Medaka model. | Per default, the workflow tries to determine the basecall model from the input data. This parameter can be used to override the detected value (or to provide a model name if none was found in the inputs). However, users should only do this if they know for certain which model was used as selecting the wrong option might give sub-optimal results. A list of recent models can be found here: https://github.com/nanoporetech/dorado#DNA-models. |  |
| medaka_model_path | string | A custom model file (.tar.gz or .hdf) to be used instead of the automatic model selection and take precedence over the optional `--override_basecaller_cfg` parameter. | Allows for users to test experimental Medaka models. Users should not provide a model with this parameter for general analysis. |  |
| large_construct | boolean | Enable assembly of larger constructs including Bacterial Artificial Chromosomes (50,000-300,000 base pairs). | Selecting this will skip approximate size filtering steps allowing the assembly of larger genomes. Multiple sequence alignment of inserts will be skipped in this mode. | False |
| trim_length | integer | Number of base pairs to trim from both ends of each read. Set to 0 if no trimming is required. |  | 0 |
| min_quality | integer | Set the minimum average quality score required for reads to be used in assembly. | Increasing this value can give more confidence in the assembly when high-quality reads are available. | 9 |
| flye_quality | string | The Flye parameter for quality of input reads, default `nano-hq`: high-quality reads, Guppy5+ SUP or Q20 (<5% error). | Other options include `nano-corr`: reads that were corrected with other methods (<3% error), `nano-raw`: pre-Guppy5 (<20% error). | nano-hq |
| non_uniform_coverage | boolean | Set this to true if your reads have highly non-uniform coverage. | Run `flye` in metagenome assembly mode, which may help with the assembly if you have high non-uniform coverage reads; generally, should not be required. | False |
| db_directory | string | Optional directory containing a gene annotation database. | A default generic annotation is provided in tar.gz format, containing entries from [fpbase](https://www.fpbase.org/), [Swiss-Prot](https://www.expasy.org/resources/uniprotkb-swiss-prot) , [Rfam](https://rfam.org/) and [snapgene](https://www.snapgene.com/) |  |
| assembly_tool | string | Select the assembly tool to use, either Canu or Flye. | Flye is the default assembler tool which will work in most cases. Alternatively select Canu but it will not work with ARM processors. | flye |
| canu_fast | boolean | Fast option can make the assembly step significantly faster. It can be used on any genome size but may produce less continuous assemblies on genomes larger than 1 Gbp. | This option is only relevant if Canu is set as the assembly_tool parameter | False |
| cutsite_mismatch | integer | Maximum number of mismatches allowed when searching for the cutsite in the reference fasta provided. Set to 0 for perfect matches only. Increasing allowed mismatches when increase risk of multiple matches, which will fail the workflow. |  | 1 |
| expected_coverage | number | The minimum coverage expected (as a percentage %) between the aligned assemblies and references if provided. This applies to both reference and assembly coverage. Applies to both the full plasmid and/or the insert. | This is used with the `--expected_identity` parameter to indicate if the construct is as expected, which is shown by a tick or cross symbol in the sample status table of the report. | 95 |
| expected_identity | number | The minimum identity expected (as a percentage %) between the aligned assemblies and references if provided. Applies to both the full plasmid and/or the insert. | This is used with the `--expected_coverage` parameter to indicate if the construct is as expected, which is shown by a tick or cross symbol in the sample status table of the report. | 99 |


### Miscellaneous Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| threads | integer | Maximum number of CPU threads to use per workflow task. | Several tasks in this workflow benefit from using multiple CPU threads. This option sets the number of CPU threads for all such processes. The total CPU resource used by the workflow is contrained by the executor configuration. | 4 |


