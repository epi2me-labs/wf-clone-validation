<!---Example of input directory structure, delete and edit as appropriate per workflow.--->
This workflow accepts either FASTQ or BAM files as input.

The FASTQ or BAM input parameters for this workflow accept one of three cases: (i) the path to a single FASTQ or BAM file; (ii) the path to a top-level directory containing FASTQ or BAM files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ or BAM files. In the first and second cases (i and ii), a sample name can be supplied with `--sample`. In the last case (iii), the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`.

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

When using a sample sheet the workflow allows the use of additional columns `approx_size`, `full_reference` `insert_reference`, `host_reference` and `regions_bedfile` which replace parameters `--approx_size`, `--full_reference`, `--insert_reference`, `--host_reference` and `--regions_bedfile` respectively. This allows per-sample variables to be applied rather than global settings. Users should provide the full path to these files, with windows users requiring to add the prefix `/mnt/c` to all paths. An example sample sheet is shown below. 

```
alias,barcode,type,approx_size,full_reference,insert_reference,host_reference,regions_bedfile
sample1,barcode01,test_sample,4000,/path/to/full_reference.fasta,/path/to/insert_reference.fasta,/path/to/host_reference.fasta,/path/to/regions_bedfile.bed
sample2,barcode02,test_sample,4000,/path/to/full_reference.fasta,/path/to/insert_reference.fasta,/path/to/host_reference.fasta,/path/to/regions_bedfile.bed
sample3,barcode03,test_sample,7000,/path/to/full_reference_alt.fasta,/path/to/insert_reference_alt.fasta/,path/to/host_reference_alt.fasta,/path/to/regions_bedfile_alt.bed
```