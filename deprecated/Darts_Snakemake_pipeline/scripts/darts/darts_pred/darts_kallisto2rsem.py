# convert kallisto tpm to rsem
# Zijun Zhang
# 7.3.2017
# revised 11.12.2017: changed to snakemake

import sys
import os
from collections import defaultdict

def read_rsem(fn):
	g2t = {}
	t2g = {}
	rsem_dict = {}
	firstline=True
	with open(fn, 'r') as fin:
		for line in fin:
			ele = line.rstrip().split()
			if firstline:
				header = {ele[x]:x for x in range(len(ele))}
				firstline = False
				continue
			gid = ele[header['gene_id']].split('.')[0]
			if not gid.startswith('ENSG'):
				continue
			tids = ele[header['transcript_id(s)']].split(',')
			tids = [x.split('.')[0] for x in tids]
			#rsem_dict[gid] = float(ele[header['TPM']])
			g2t[gid] = tids
			for tid in tids:
				t2g[tid] = gid
	return rsem_dict, t2g


def read_kallisto(fn, t2g):
	kal_dict = defaultdict(float)
	with open(fn, 'r') as fin:
		firstline=True
		for line in fin:
			ele = line.rstrip().split()
			if firstline:
				header = {ele[x]:x for x in range(len(ele))}
				firstline=False
				continue
			tid = ele[header['target_id']].split('.')[0]
			tpm = float(ele[header['tpm']])
			kal_dict[t2g[tid]] += tpm
	return kal_dict
	
def write_out_comparison(rsem_dict, kal_dict, fn):
	with open(fn, 'w') as fout:
		fout.write('gene_id\trsem\tkallisto\n')
		for g in kal_dict:
			fout.write('%s\t%f\t%f\n'%(g, rsem_dict[g], kal_dict[g]))
	return


def write_out_quant(kal_dict, fn):
	with open(fn, 'w') as fout:
		fout.write('gene_id\tTPM\n')
		for g in kal_dict:
			fout.write('%s\t%f\n'%(g, kal_dict[g]))
	return

def kallisto2rsem(kal_dir, t2g):
	kal_fn = os.path.join(kal_dir, 'abundance.tsv')
	kal_dict = read_kallisto(kal_fn, t2g)
	kal_outfn = os.path.join(kal_dir, 'gene_tpm.tsv')
	write_out_quant(kal_dict, kal_outfn)
	return


def parser(kal_dir):
	srcdir = os.path.dirname(os.path.realpath(__file__))
	#rsem_fn = '/u/nobackup/yxing/NOBACKUP/frankwoe/ENCODE/CRISPR_RNAseq_201707/PTBP1_K562_ENCSR415DJT/rsem/genequant_rep1_ENCFF589WON.tsv'
	rsem_fn = os.path.join(srcdir, 'rsem_template.tsv')
	rsem_dict, t2g = read_rsem(rsem_fn)
	kallisto2rsem(kal_dir, t2g)


