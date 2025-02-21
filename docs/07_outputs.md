Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| worfklow report | wf-clone-validation-report.html | A report bringing together the main results of the workflow, across samples. | aggregated |
| sample status | sample_status.txt | A CSV file with per-sample assembly success or failure reasons | aggregated |
| plasmid annotations | plannotate.json | Plasmid annotations in a JSON structure. | aggregated |
| annotations bed | {{ alias }}.annotations.bed | Plasmid annotations in a BED file format for onward use | per-sample |
| annotations gbk | {{ alias }}.annotations.gbk | Plasmid annotations in a GBK file format for onward use | per-sample |
| Assembly FASTQ | {{ alias }}.final.fastq | Sequence and quality score of the final assembly. | per-sample |
| Assembly statistics | {{ alias }}.assembly_stats.tsv | Assembly statistics from fastcat. | per-sample |
| Insert FASTA | {{ alias }}.insert.fasta | Insert sequence found in the final assembly, only relevant if the primers parameter was used. | per-sample |
| Variant stats report | {{ alias }}.full_construct.stats | A BCF stats report with any variants found, only relevant if a full reference was provided. | per-sample |
| Variants BCF file | {{ alias }}.full_construct.calls.bcf | A BCF file with any variants found per sample, only relevant if a full reference was provided. | per-sample |
| Reference alignment | {{ alias }}.bam | Reference aligned with the assembly in BAM format, only relevant if a full reference was provided. | per-sample |
| Reference alignment index | {{ alias }}.bam.bai | The index for the reference aligned with the assembly, only relevant if a full reference was provided. | per-sample |
| Host reference alignment | {{ alias }}.host.bam | Host reference aligned with sample in BAM format, only relevant if a host reference was provided. | per-sample |
| Host reference alignment index | {{ alias }}.host.bam.bai | The index for the host reference aligned with sample, only relevant if a host reference was provided. | per-sample |
| BAM Stats | {{ alias }}.bam.stats | Stats report for the reference aligned with the assembly, only relevant if a full reference was provided. | per-sample |
