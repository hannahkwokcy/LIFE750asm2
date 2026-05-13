# LIFE750asm2
This repository contains the codes, intermediate files and results of the LIFE750 assignment2.
# Data acquisition
The dataset was obtained via the linx command: wget --mirror --no-parent -nd -R "index*" https://cgr.liv.ac.uk/454/acdarby/LIFE750/Life750_datasets/dataset_hlckwok2/
# Description
* `1_variant_calling`: Contains the resulting .vcf files of variant calling and the commands used for this part.
* `2_transcriptomics`: Contains the R script for transcriptomic analysis and graphs generation. This also contains the summary differential expression results file as well as heatmap, PCA and volcano plot graphs generated from the code.
* `3_genome_binding`: Contains the command used and resulting gene binding summary file.
* `gene_x`: Contains the orignal files retrieved during data acquisition as well as intermediate files generated during variant calling.
* `annotations`, `cutrun`, `expression`, `gene_x`, `metadata`: Contains only the original files retrieved during the data acquisition step.
