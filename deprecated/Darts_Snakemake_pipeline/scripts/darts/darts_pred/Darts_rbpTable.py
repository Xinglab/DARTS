# read in folder structure, prepare the data format,
# and call Darts[Spider] to run on processed splicing data
# Zijun Zhang
# 5.18.2017
# revised 7.1.2017: use only one replicate as downsampling
# revised 7.2.2017: enable replicate-wise or merged one-click pipeline
# revised 7.3.2017: remove replicate-wise, add random down-sampling
# revised 8.30.2017: changed to Roadmap data
# revised 11.12.2017: changed to snakemake package

import sys
import os
from collections import defaultdict
from datetime import datetime
import numpy as np


def read_kal_fn(tsv_fn_list):
	gene_exp = defaultdict(list)
	for genequant in tsv_fn_list:
		#print(genequant)
		with open(genequant,'r') as f:
			firstline=True
			for line in f:
				ele = line.rstrip().split()
				if firstline:
					header = {ele[x]:x for x in range(len(ele))}
					firstline=False
					continue
				gene_id = ele[header['gene_id']].split('.')[0]
				#if gene_id in RBP:
				gene_exp[gene_id].append(float(ele[header['TPM']]))
	for gene in gene_exp:
		gene_exp[gene] = np.mean(gene_exp[gene])
	return gene_exp


def upper_quantile_normalization(gene_exp):
	tpm = np.asarray(gene_exp.values())
	uq = np.percentile(tpm[tpm!=0], 75)
	for gene in gene_exp:
		gene_exp[gene] /= uq
	return gene_exp


def read_rbp_list(rbp_fn):
	RBP = []
	RBP_geneName = {}
	with open(rbp_fn, 'r') as f:
		for line in f:
			ele = line.rstrip().split()
			RBP.append(ele[0])
			RBP_geneName[ele[0]] = ele[1]
	return RBP, RBP_geneName	


	
def make_rbp_exp_table(target, control, darts_dir, RBP, RBP_geneName, s1, s2):
	#global beta0, beta1
	tar_gene_exp = read_kal_fn(target)
	#tar_gene_exp = upper_quantile_normalization(tar_gene_exp)
	con_gene_exp = read_kal_fn(control)
	#con_gene_exp = upper_quantile_normalization(con_gene_exp)
	rbp_exp_fn = os.path.join(darts_dir, 'RBP_tpm.txt')
	with open(rbp_exp_fn, 'w') as fout:
		fout.write('\t%s\t%s\n'%(s1, s2))
		for rbp in RBP:
			if rbp in tar_gene_exp:
				tar_tpm = tar_gene_exp[rbp]
				#tar_tpm = tar_tpm*beta1+beta0
			else:
				tar_tpm = 'NA'
				print("NA produced")
			if rbp in con_gene_exp:
				con_tpm = con_gene_exp[rbp]
				#con_tpm = con_tpm*beta1+beta0
				#con_tpm = 0 if con_tpm<0 else con_tpm
			else:
				con_tpm = 'NA'
				print("NA produced")
			line = '{0}\t{1}\t{2}\n'.format(RBP_geneName[rbp], tar_tpm, con_tpm)
			fout.write(line)
	return


def parser(outdir, comparison, kallisto_1, kallisto_2):
	srcdir = os.path.dirname(os.path.realpath(__file__))
	rbp_fn = os.path.join(srcdir, 'Tuschl_RBP_list.kallisto.txt')
	RBP, RBP_geneName = read_rbp_list(rbp_fn)
	s1, s2 = comparison.split('-')
	darts_dir = outdir
	make_rbp_exp_table(kallisto_1, kallisto_2, darts_dir, RBP, RBP_geneName, s1, s2)


if __name__ == '__main__':
	write_dir = '/u/nobackup/yxing/NOBACKUP/frankwoe/Roadmap/Darts_100_101'
	sample_fn = '/u/nobackup/yxing/NOBACKUP/frankwoe/Roadmap/sample.readLen_100_101.txt'
	rbp_fn = '/u/nobackup/yxing/NOBACKUP/frankwoe/Roadmap/Darts_caller/Tuschl_RBP_list.kallisto.txt'
	sample_dict = read_sample_info(sample_fn)
	RBP, RBP_geneName = read_rbp_list(rbp_fn)
	make_rbpexp_table_exhaustive(sample_dict, write_dir, RBP, RBP_geneName)	
