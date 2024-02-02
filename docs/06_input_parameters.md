### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| fastq | string | FASTQ files to use in the analysis. | This accepts one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| approx_size | integer | Approximate size of the plasmid in base pairs. This can also be defined per sample: see sample_sheet parameter. |  | 7000 |
| assm_coverage | integer | Fold coverage for use per assembly | This is the coverage that will be used to subsample reads to use for the assembly. | 60 |
| primers | string | TSV File containing primers used to find inserts. If left empty then inserts will not be searched for. | Specify one or more primer sets which will be used to find the sequence inserted in the construct. This file should be in .tsv format containing columns [primer_name,  5' primer, 3' primer] with no header. An example `primers.tsv` for pRham/T7 is available in the demo data for the workflow. |  |
| analyse_unclassified | boolean | Analyse unclassified reads from input directory. By default the workflow will not process reads in the unclassified directory. | If selected and if the input is a multiplex directory the workflow will also process the unclassified directory. | False |
| basecaller_cfg | string | Name of the model that was used to basecall signal data, used to select an appropriate Medaka model. | The basecaller configuration is used to automatically select the appropriate Medaka model. The automatic selection can be overridden with the 'medaka_model' parameter. The model list only shows models that are compatible with this workflow. | dna_r10.4.1_e8.2_400bps_sup@v4.2.0 |


### Reference Genome Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| insert_reference | string | Optional file containing insert reference sequence which will be used for comparison with consensus insert in the report. | Providing a reference sequence can be useful as a QC on the base-level resolution of the the reconstructed insert sequences. |  |
| host_reference | string | A host reference genome FASTA file. Read which map to this reference are discarded and not used for the assembly. |  |  |
| regions_bedfile | string | If a host_reference supplied, add an optional BED file to provide host reference regions that will be masked during filtering. |  |  |


### Sample Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_sheet | string | A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a directory containing sub-directories with FASTQ files. An optional column `approx_size` can be added to provide size estimates for each sample. When not provided, the `--approx_size` parameter will be used for all samples. | The sample sheet is a CSV file with, minimally, columns named `barcode` and `alias`. Extra columns are allowed. A `type` column is required for certain workflows and should have the following values; `test_sample`, `positive_control`, `negative_control`, `no_template_control`. |  |
| sample | string | A single sample name for non-multiplexed data. Permissible if passing a single .fastq(.gz) file or directory of .fastq(.gz) files. |  |  |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all user-facing files. |  | output |
| prefix | string | The prefix attached to each of the output filenames. |  |  |


### Advanced Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| large_construct | boolean | Enable assembly of larger constructs including Bacterial Artificial Chromosomes (50,000-300,000 base pairs). | Selecting this will skip approximate size filtering steps allowing the assembly of larger genomes. Multiple sequence alignment of inserts will be skipped in this mode. | False |
| trim_length | integer | Number of base pairs to trim from the both ends of read. |  | 150 |
| medaka_model | string | The name of a Medaka model to use. By default the workflow will select an appropriate Medaka model from the basecaller configuration provided. Entering a name here will override the automated selection and use the Medaka model named here. | The workflow will attempt to map the basecalling model used to a suitable Medaka model. You can override this by providing a model with this option instead. |  |
| flye_quality | string | The Flye parameter for quality of input reads, default `nano-hq`: high-quality reads, Guppy5+ SUP or Q20 (<5% error). | Other options include `nano-corr`: reads that were corrected with other methods (<3% error), `nano-raw`: pre-Guppy5 (<20% error). | nano-hq |
| non_uniform_coverage | boolean | Set this to true if your reads have highly non-uniform coverage. | Run `flye` in metagenome assembly mode, which may help with the assembly if you have high non-uniform coverage reads; generally, should not be required. | False |
| db_directory | string | Optional directory containing a gene annotation database. | A default generic annotation is provided in tar.gz format, containing entries from [fpbase](https://www.fpbase.org/), [Swiss-Prot](https://www.expasy.org/resources/uniprotkb-swiss-prot) , [Rfam](https://rfam.org/) and [snapgene](https://www.snapgene.com/) |  |


### Miscellaneous Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| threads | integer | Maximum number of CPU threads to use per workflow task. | Several tasks in this workflow benefit from using multiple CPU threads. This option sets the number of CPU threads for all such processes. The total CPU resource used by the workflow is contrained by the executor configuration. | 4 |
| disable_ping | boolean | Enable to prevent sending a workflow ping. |  | False |


