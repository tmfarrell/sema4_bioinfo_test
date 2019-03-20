#!/usr/bin/env python3 

##
## retrieve_coding_exons_for_genes.py
## 
## Tim Farrell (tfarrell01@gmail.com)
## 20190319
## 

import argparse
import pandas as pd
import mysql.connector 

## FUNCTIONS 
# intervals overlap? 
def do_overlap(i1_start, i1_end, i2_start, i2_end):
	if i2_start <= i1_start and i1_start <= i2_end: 
		return(True)
	if i2_start <= i1_end and i1_end <= i2_end: 
		return(True)
	return(False)


## MAIN 
def main(): 
	# parse cmd-line args
	parser = argparse.ArgumentParser()
	parser.add_argument('--genes', nargs='+', 
		help='List of gene names to retrieve coding exons for.')
	args = parser.parse_args()

	# open connection to db, and init cursor 
	cnx = mysql.connector.connect(user='genome', database='hg38', host='genome-mysql.cse.ucsc.edu', port=3306)
	crs = cnx.cursor(buffered=True)

	gene_exons_df = pd.DataFrame()
	# for each gene 
	for gene in args.genes: 
		# get best ncbiRefSeq id, with valid protein accession (i.e. coding)
		#print('querying genome-mysql.cse.ucsc.edu for gene: %s...' % gene)
		crs.execute('SELECT id FROM ncbiRefSeqLink WHERE name = "%s" and protAcc != ""' % gene)
		ncbi_id = crs.fetchone()[0]
		# retrieve the predicted exons for that id 
		crs.execute('SELECT transcriptId, seqId, exonId, chrom, chromStart, chromEnd FROM wgEncodeGencodeExonSupportV28 WHERE seqId = "%s"' % ncbi_id)
		# init new dataframe for gene
		df = pd.DataFrame(crs.fetchall(), columns=['transcipt_id','ncbi_id','exon_id','chrom','start','end'])
		# drop coordinate duplicates, and sort 
		#print('resolving overlapping exons for gene: %s...' % gene)
		df = (df.loc[df[['chrom','start','end']].drop_duplicates().index, :].sort_values(['chrom','start','end']).reset_index(drop=True))
		# check that no intervals overlap 
		non_overlapping_index = []
		for i in range(len(df.index) - 1):
			idx_i = df.index[i]
			idx_j = df.index[i+1]
			chrom_i, start_i, end_i = df.loc[idx_i, ['chrom','start','end']]
			chrom_j, start_j, end_j = df.loc[idx_j, ['chrom','start','end']]
			if chrom_i != chrom_j or not do_overlap(start_i, end_i, start_j, end_j):
				non_overlapping_index = [i for i in non_overlapping_index if i != idx_i and i != idx_j] + [idx_i, idx_j]
			else:
				non_overlapping_index = [i for i in non_overlapping_index if i != idx_i] + [idx_j]
		# add to genes dataframe 
		gene_exons_df = gene_exons_df.append(df.loc[non_overlapping_index,:], ignore_index=True)

	# print coordinates in bed format to STDOUT
	print(gene_exons_df[['chrom','start','end']].to_csv(sep='\t', index=False, header=False))
	return  



if __name__ == '__main__': 
	main() 