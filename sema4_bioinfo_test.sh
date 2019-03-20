#!/usr/bin/env bash 

pip3 install -r requirements.txt

# question 1
python3 compare_two_fastas.py --help 
python3 compare_two_fastas.py --fasta_refs \
	http://hgdownload.cse.ucsc.edu/goldenPath/hg19/hg19Patch13/hg19Patch13.fa.gz \
	http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz


# question 2 
python3 retrieve_coding_exons_for_genes.py --help
python3 retrieve_coding_exons_for_genes.py --genes BRCA1 BRCA2


# question 3 
Rscript vcf_to_mutation_types_table.R --help 
Rscript vcf_to_mutation_types_table.R --vcf Q3.vcf --mutation_types_table mutation_types.txt 