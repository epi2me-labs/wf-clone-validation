## Introduction

The workflow accepts FASTQ as the primary input.

Optional inputs:

* A `host_reference` FASTA file.
  
* A primers TSV file if they differ from the default pRham/T7.

The steps of the workflow are as follows - 

* If a `host_reference` is provided, [minimap2](https://github.com/lh3/minimap2) is used to align the reads to the `host_reference`, aligned reads are filtered out.

* The reads are then trimmed at the ends using [Seqkit](https://bioinf.shenwei.me/seqkit/) with the provided `--trim_length` parameter. 

* The sequences are then downsampled using the tool [Rasusa](https://github.com/mbhall88/rasusa).

* [Trycycler](https://github.com/rrwick/Trycycler) is used to create 3 subsamples which are each assembled using [Flye](https://github.com/fenderglass/Flye).

* If there are concatemers in the assembly these are found using minimap2 and deconcatenated using a Python script. 

* Trycycler is used to reconcile the subsampled assemblies into one final assembly. This is polished with [Medaka](https://github.com/nanoporetech/medaka).

* [Seqkit](https://bioinf.shenwei.me/seqkit/) is used to find inserts using the primers supplied to the `--primers` parameter. 

* The assembly is annotated using [pLannotate](https://github.com/barricklab/pLannotate) with the default database containing entries from [fpbase](https://www.fpbase.org/), [Swiss-Prot](https://www.expasy.org/resources/uniprotkb-swiss-prot), [Rfam](https://rfam.org/) and [snapgene](https://www.snapgene.com/). 

* A quality score for the assembly is provided by Medaka.

* Optionally a reference insert sequence can be provided which is aligned to the consensus and any variants are reported by [bcftools](https://samtools.github.io/bcftools/bcftools.html).