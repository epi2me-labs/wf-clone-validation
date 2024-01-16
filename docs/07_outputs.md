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
