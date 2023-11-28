*I don't know the approximate size of my plasmid?* - On most occasions you can use the mode of the data as an approximate guide To find the mode, you can run the workflow with the default settings, and from the raw data read length plot find the highest peak. This value should approximate the plasmid size because for most plasmids only one cut is made to the circular plasmid prior to sequencing, meaning each read is of the full plasmid. Furthermore, it is better to overestimate the approximate size than underestimate.

*Does the workflow report contaminants?* - The workflow has no way of reporting contaminants. However, if contaminants are present, the workflow may struggle to create consistent assemblies and the output assemblies are likely to show low quality. If you have a reference for an expected contaminant, you could use this as the host reference to filter out any reads that align with that.

*Can I use my own annotation database?* – Currently using your own annotation database is not supported, but we may add it in future.

*Does this workflow support reference based assembly?* - It does not have a reference based assembly mode.

*Does this workflow have support for bacterial artificial chromosomes (BACs)?* - This workflow does not yet have BAC support and has not been tested for assembly of genomes larger than 50,000bps

If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-template/issues) page or start a discussion on the [community](https://community.nanoporetech.com/).