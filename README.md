# Clone validation workflow

*De novo* assembly of plasmid sequences, designed for verifying the results of molecular cloning experiments.



## Introduction

Among other uses this workflow could determine the success of a molecular cloning experiment and determine whether one DNA sequence has been correctly inserted into another as an experimentalist was expecting.

In brief, this workflow will perform the following:

+ *De novo* assembly of plasmids.
+ Annotation of the full assembly.
+ Provide a per base quality score of the plasmid assembly.
+ Locate an insert sequence in a plasmid using provided primers.
+ Multiple sequence alignment between insert sequences from different samples.
+ Create an assembly dot plot showing repetitive regions in the created assemblies. 
+ Comparison between an insert reference and the assembled insert.



## Compute requirements

Recommended requirements:

+ CPUs = 4
+ Memory = 8GB

Minimum requirements:

+ CPUs = 4
+ Memory = 8GB

Approximate run time: 6 minutes per sample for 10,000 reads

ARM processor support: True




## Install and run

These are instructions to install and run the workflow on command line. You can also access the workflow via the [EPI2ME application](https://labs.epi2me.io/downloads/).

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and software resources, therefore nextflow will need to be installed before attempting to run the workflow.

The workflow can currently be run using either [Docker](https://www.docker.com/products/docker-desktop) or
[singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker or singularity is installed. This is controlled by the [`-profile`](https://www.nextflow.io/docs/latest/config.html#config-profiles) parameter as exemplified below.

It is not required to clone or download the git repository in order to run the workflow.
More information on running EPI2ME workflows can be found on our [website](https://labs.epi2me.io/wfindex).

The following command can be used to obtain the workflow. This will pull the repository in to the assets folder of nextflow and provide a list of all parameters available for the workflow as well as an example command:

```
nextflow run epi2me-labs/wf-clone-validation –help
```
A demo dataset is provided for testing of the workflow. It can be downloaded using:
```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-clone-validation/clone_val_test.tar.gz
tar -xzvf clone_val_test.tar.gz
```
The workflow can be run with the demo data using:
```
nextflow run epi2me-labs/wf-clone-validation \
--fastq clone_val_test/fastq --primers clone_val_test/primers.tsv \
--host_reference clone_val_test/host_reference.fa.gz --regions_bedfile clone_val_test/reference.bed \
--insert_reference clone_val_test/insert_reference.fasta --sample_sheet clone_val_test/sample_sheet.csv \
-profile standard
```
For further information about running a workflow on the cmd line see https://labs.epi2me.io/wfquickstart/




## Related protocols

This workflow is designed to take input sequences that have been produced from [Oxford Nanopore Technologies](https://nanoporetech.com/) devices using this protocol:

https://community.nanoporetech.com/docs/prepare/library_prep_protocols/plasmid-sequencing-using-sqk-rbk004/



## Input example

<!---Example of input directory structure, delete and edit as appropriate per workflow.--->
This workflow accepts FASTQ files as input.

The FASTQ input parameter for this workflow accepts one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second cases (i and ii), a sample name can be supplied with `--sample`. In the last case (iii), the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`.

```
(i)                     (ii)                 (iii)    
input_reads.fastq   ─── input_directory  ─── input_directory
                        ├── reads0.fastq     ├── barcode01
                        └── reads1.fastq     │   ├── reads0.fastq
                                             │   └── reads1.fastq
                                             ├── barcode02
                                             │   ├── reads0.fastq
                                             │   ├── reads1.fastq
                                             │   └── reads2.fastq
                                             └── barcode03
                                              └── reads0.fastq
```



## Input parameters

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
| trim_length | integer | Number of base pairs to trim from both ends of each read. |  | 150 |
| medaka_model | string | The name of a Medaka model to use. By default the workflow will select an appropriate Medaka model from the basecaller configuration provided. Entering a name here will override the automated selection and use the Medaka model named here. | The workflow will attempt to map the basecalling model used to a suitable Medaka model. You can override this by providing a model with this option instead. |  |
| flye_quality | string | The Flye parameter for quality of input reads, default `nano-hq`: high-quality reads, Guppy5+ SUP or Q20 (<5% error). | Other options include `nano-corr`: reads that were corrected with other methods (<3% error), `nano-raw`: pre-Guppy5 (<20% error). | nano-hq |
| non_uniform_coverage | boolean | Set this to true if your reads have highly non-uniform coverage. | Run `flye` in metagenome assembly mode, which may help with the assembly if you have high non-uniform coverage reads; generally, should not be required. | False |
| db_directory | string | Optional directory containing a gene annotation database. | A default generic annotation is provided in tar.gz format, containing entries from [fpbase](https://www.fpbase.org/), [Swiss-Prot](https://www.expasy.org/resources/uniprotkb-swiss-prot) , [Rfam](https://rfam.org/) and [snapgene](https://www.snapgene.com/) |  |
| assembly_tool | string | Select the assembly tool to use, either Canu or Flye. | Flye is the default assembler tool which will work in most cases. Alternatively select Canu but it will not work with ARM processors. | flye |
| canu_fast | boolean | Fast option can make the assembly step significantly faster. It can be used on any genome size but may produce less continuous assemblies on genomes larger than 1 Gbp. | This option is only relevant if Canu is set as the assembly_tool parameter | False |


### Miscellaneous Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| threads | integer | Maximum number of CPU threads to use per workflow task. | Several tasks in this workflow benefit from using multiple CPU threads. This option sets the number of CPU threads for all such processes. The total CPU resource used by the workflow is contrained by the executor configuration. | 4 |
| disable_ping | boolean | Enable to prevent sending a workflow ping. |  | False |






## Outputs

Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| worfklow report | ./wf-clone-validation-report.html | A report bringing together the main results of the workflow, across samples. | aggregated |
| sample status | ./sample_status.txt | A CSV file with per-sample assembly success or failure reasons | aggregated |
| plasmid annotations | ./plannotate.json | Plasmid annotations in a JSON structure. | aggregated |
| annotations bed | ./{{ alias }}.annotations.bed | Plasmid annotations in a BED file format for onward use | per-sample |
| annotations gbk | ./{{ alias }}.annotations.gbk | Plasmid annotations in a GBK file format for onward use | per-sample |
| Assembly FASTQ | ./{{ alias }}.final.fastq | Sequence and quality score of the final assembly. | per-sample |
| Insert FASTA | ./{{ alias }}.insert.fasta | Insert sequence found in the final assembly, only relevant if the primers parameter was used. | per-sample |




## Pipeline overview

### 1. Concatenates input files and generate per read stats.

The [fastcat](https://github.com/epi2me-labs/fastcat) tool is used to concatenate multifile samples to be processed by the workflow. It will also output per read stats including average read lengths and qualities. Reads with lengths between 0.5 X and 1.5 X the approximate size are filtered out unless the `--large_construct` parameter is provided which indicates the assembly is expected to be larger (50,000-300,000 bps).

### 2. Filter out host reference reads

If a host_reference fasta file is provided, [Minimap2](https://github.com/lh3/minimap2) is used to align all reads to the host_reference, and any aligned reads are filtered out.

### 3. Trim reads

The reads are then trimmed at the ends using [SeqKit](https://bioinf.shenwei.me/seqkit/) with the provided trim length parameter, which has a default of 150bps. This default was chosen as we found it to produce optimum results and we don’t recommend changing this. At this stage SeqKit is also used to filter out reads that are greater than 1.2 x the approximate size or less than 100bps.

### 4. Subsample reads

The sequences are then subsampled using [Rasusa](https://github.com/mbhall88/rasusa). The subsampling will take the expected coverage parameter in to account; as we will be repeating the assembly 3 times, we subsample to a target of approximately 3x the required coverage. However, this is just a target and if there is not enough data, Rasusa is still able to create the 3 subsamples. The approximate size parameter is also used by Rasusa to work out the target number of bases and therefore number of reads required for each of the subsamples.

### 5. Create 3 subsamples

+[Trycycler](https://github.com/rrwick/Trycycler) is used to create 3 subsamples as we will be creating three assemblies and finding the consensus between all three. This consensus generation will be handled by Ttrycycler.

### 6. Assembly

We perform the assembly for each of the 3 subsamples separately. The assembly is done using either [Flye](https://github.com/fenderglass/Flye) or [Canu](https://github.com/marbl/canu) depending on what is set as the `assembly_tool` parameter. Both Flye and Canu are popular assemblers that usually produce reliable assemblies. Flye is our default assembler as it usually provides reliable assemblies in less time than Canu, and supports ARM processors. If Flye fails to assemble you may wish to try Canu.

### 7. De-concatenate

If there are concatemers in the assembly, these are found using minimap2 and de-concatenated using a custom Python script. If the assembly is already roughly the expected approximate size, this de-concatenate step will be skipped.

### 8. Reconcile and polish

Trycycler is used to reconcile the subsampled assemblies into one final assembly. This is then polished with [Medaka](https://github.com/nanoporetech/medaka). A per-base quality score for the assembly is output by Medaka in a FASTQ file. This is used for creating the mean assembly quality you will find in the report.

### 8. Insert location and QC

SeqKit is used to locate inserts using the primers supplied to the primers parameter.

A multiple sequence alignment (MSA) will be done using [Pyspoa](https://github.com/nanoporetech/pyspoa). This will be presented in the report to help users compare inserts across samples in a multi-sample run. If an insert reference FASTA file is provided, this will also be included in the MSA.

If a reference insert FASTA sequence is provided, [BCFtools](https://samtools.github.io/bcftools/bcftools.html) is used to find variants between it and the final insert assembly, and are reported in BCF file per sample.

### 10. Annotate

The assembly is annotated by pLannotate](https://github.com/barricklab/pLannotate) to show any features that are present. The default database is used, which contains entries from [FPbase](https://www.fpbase.org/), [Swiss-Prot](https://www.expasy.org/resources/uniprotkb-swiss-prot), [Rfam](https://rfam.org/) and [SnapGene](https://www.snapgene.com/). Descriptions, percentage match and length of the match are also provided.

### 11. Self alignment

For each sample a self alignment will be done using [Last](https://gitlab.com/mcfrith/last) and the output will be presented as a dotplot. This can help identify any repetitive regions in your final assembly.



## Troubleshooting

+ If there are no assemblies output by the workflow, open the wf-clone-validation-report.html to look at failure reasons. Check the read summary section for quality and quantity
of reads before and after downsampling to ensure there is enough data for the assembly. If there is not sufficient data, you may need to adjust the approx size and coverage options.
+ If the workflow fails please run it with the demo data set to ensure the workflow itself is working. This will help us determine if the issue is related to the environment, input parameters or a bug.
+ See how to interpret some common nextflow exit codes [here](https://labs.epi2me.io/trouble-shooting/).



## FAQ's

*I don't know the approximate size of my plasmid?* - On most occasions you can use the mode of the data as an approximate guide To find the mode, you can run the workflow with the default settings, and from the raw data read length plot find the highest peak. This value should approximate the plasmid size because for most plasmids only one cut is made to the circular plasmid prior to sequencing, meaning each read is of the full plasmid. Furthermore, it is better to overestimate the approximate size than underestimate.

*Does the workflow report contaminants?* - The workflow has no way of reporting contaminants. However, if contaminants are present, the workflow may struggle to create consistent assemblies and the output assemblies are likely to show low quality. If you have a reference for an expected contaminant, you could use this as the host reference to filter out any reads that align with that.

*Can I use my own annotation database?* – Currently using your own annotation database is not supported, but we may add it in future.

*Does this workflow support reference based assembly?* - It does not have a reference based assembly mode.

*Does this workflow have support for bacterial artificial chromosomes (BACs)?* - This workflow does not yet have BAC support and has not been tested for assembly of genomes larger than 50,000bps

If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-template/issues) page or start a discussion on the [community](https://community.nanoporetech.com/).



## Related blog posts

+ [Assembly tools and Flye](https://labs.epi2me.io/assembly-flye/)
 
 See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts.



