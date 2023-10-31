#!/bin/bash
##### VELOCYTO ####
# run on each sample, refdata-gex-mm10-2020-A used as reference

samtools sort -t CB -O BAM -o cellsorted_possorted_genome_bam.bam $file_bam
velocyto run -b barcodes.tsv -o $output_path -m mm10_rmsk.gtf $file_bam genes.gtf
