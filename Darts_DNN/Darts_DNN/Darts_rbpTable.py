"""
Darts_DNN - rbpTable

Utilities for reading and preparinig trans features from gene expression
"""
import sys
import os
from collections import defaultdict
from datetime import datetime
import numpy as np
from . import config


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


def read_rbp_list(rbp_str):
	RBP = []
	RBP_geneName = {}
	#with open(rbp_fn, 'r') as f:
	if 1:
		#for line in f:
		for line in rbp_str.split('\n'):
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


def parser(outdir, kallisto_1, kallisto_2):
	rbp_str = open(config.RBP_GENE_LIST_PATH, 'r').read()
	RBP, RBP_geneName = read_rbp_list(rbp_str)
	darts_dir = outdir
	make_rbp_exp_table(kallisto_1, kallisto_2, darts_dir, RBP, RBP_geneName, 'condition_1', 'condition_2')

