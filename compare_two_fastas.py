#!/usr/bin/env python3

##
## compare_ref_fastas.py 
## 
## Tim Farrell (tfarrell01@gmail.com)
## 20190319
## 

# import libraries 
import os 
import wget 
import gzip
import shutil
import argparse 
import itertools
from Bio import SeqIO


def main(): 
	# build cmd-line parser  
	parser = argparse.ArgumentParser()
	parser.add_argument('--fasta_refs', required=True, nargs=2, 
		help="Local paths or urls to 2 reference fastas (if urls, they will be downloaded using wget).")

	# parse cmd-line args 
	args = parser.parse_args()

	# open fasta indices 
	try:
		# try to open as local files
		f_idx = [] 
		for i in range(2): 
			basename, ext = os.path.splitext(args.fasta_refs[i])
			if ext == '.gz':
				with gzip.open(args.fasta_refs[i], 'rb') as f_in:
					with open(basename, 'wb') as f_out:
						print("decompressing %s..." % args.fasta_refs[i])
						shutil.copyfileobj(f_in, f_out)
				print("reading index of %s..." % basename)
				f_idx.append(SeqIO.index(basename, "fasta"))
			else:
				print("reading index of %s..." % args.fasta_refs[i])
				f_idx.append(SeqIO.index(args.fasta_refs[i], "fasta"))
	except: 	
		# else attempt to download from url, and then open 
		f_idx = [] 
		for i in range(2): 
			print('downloading %s...' % args.fasta_refs[i])
			local_fasta_ref = './' + os.path.basename(args.fasta_refs[i])
			wget.download(args.fasta_refs[i], local_fasta_ref)
			# gunzip, if necessary, then open index
			basename, ext = os.path.splitext(local_fasta_ref)
			if ext == '.gz':
				with gzip.open(local_fasta_ref, 'rb') as f_in:
					with open(basename, 'wb') as f_out:
						print("decompressing %s..." % local_fasta_ref)
						shutil.copyfileobj(f_in, f_out)
				print("reading index of %s..." % basename)
				f_idx.append(SeqIO.index(basename, "fasta"))
			else: 
				print("reading index of %s..." % local_fasta_ref)
				f_idx.append(SeqIO.index(local_fasta_ref, "fasta"))
	
	# get lists of common and unique chromosomes from each ref fasta index 
	common_chroms = [c0 for c0, c1 in itertools.product(f_idx[0].keys(), f_idx[1].keys()) if c0 == c1]
	f0_unique_chroms = list(set(f_idx[0].keys()).difference(set(f_idx[1].keys())))
	f1_unique_chroms = list(set(f_idx[1].keys()).difference(set(f_idx[0].keys())))

	# loop through commmon chroms, testing and printing chromosome equality
	for c in common_chroms: 
		chrom0 = f_idx[0][c].seq
		chrom1 = f_idx[1][c].seq
		print(c, 'same' if chrom0 == chrom1 else 'different')
	
	# then print out those that are missing from each, respectively
	for c in f0_unique_chroms: 
		print(c, 'missing from %s' % args.fasta_refs[1])
	for c in f1_unique_chroms: 
		print(c, 'missing from %s' % args.fasta_refs[0])

	return 


if __name__ == "__main__": 
	main()

