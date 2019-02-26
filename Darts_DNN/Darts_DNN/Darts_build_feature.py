# -*- coding: UTF-8 -*-

"""
Darts_DNN - Darts_build_feature

make the sequence and RBP features for a given comparison
"""

import os
import sys
import pandas as pd
import numpy as np
import h5py
from collections import defaultdict
import logging
logger = logging.getLogger('Darts_DNN.build_feature')

from . import Darts_rbpTable as _Darts_rbp
from . import Darts_tx2g as _Darts_tx2g
from . import config
from .utils import read_sequence_feature

def read_sp_out(sp_fn, sig_only=False):
	if not os.path.isfile(sp_fn):
		raise Exception("Darts-flat result not found: %s"%sp_fn)
		return None
	data=pd.read_table(sp_fn)
	data.index = data.ID
	if sig_only:
		evt = np.asarray(data.ID[[i for i in range(data.shape[0]) if data.post_pr[i]>0.9 or data.post_pr[i]<0.1]])
		label = np.asarray([1 if data.loc[x,'post_pr']>0.9 else 0 for x in evt])
	else:
		evt = np.asarray(data.ID[[i for i in data.ID
			if not np.isnan(data.post_pr[i])]])
		label = np.asarray([data.loc[x,'post_pr'] for x in evt])
	metadata = np.array([[x, 'condition_1', 'condition_2', '', data.I1[x], data.S1[x], data.I2[x], data.S2[x]] for x in evt], dtype=str).astype('|S60')
	return {'evt':evt, 'label':label, 'metadata':metadata}


def read_rbp_exp(rbp_fn):
	colnames = np.empty((config.TRANS_TOTAL_NUM*2), dtype='str').astype('|S30')
	rbp_exp_vec = np.empty((config.TRANS_TOTAL_NUM*2))
	with open(rbp_fn, 'r') as fin:
		idx = 0
		firstline = True
		for line in fin:
			if firstline:
				firstline=False
				continue
			ele = line.rstrip().split()
			colnames[idx] = 'KD_' + ele[0]
			colnames[idx+config.TRANS_TOTAL_NUM] = 'Con_' + ele[0]
			rbp_exp_vec[idx] = float(ele[1])
			rbp_exp_vec[idx+config.TRANS_TOTAL_NUM] = float(ele[2])
			idx += 1
	return colnames, rbp_exp_vec


def read_current_rbp_max(fn):
	return open(fn, 'r').read()


def RBPexp_normalizer(rbp_exp_vec, event_type):
	RBPexp_absmax = np.zeros((config.TRANS_TOTAL_NUM))
	rbp_max_str = read_current_rbp_max(config.CURRENT_TRANS_PATH[event_type])
	if 1:
		i=0
		for line in rbp_max_str.split('\n'):
			try:
				RBPexp_absmax[i] = float(line.rstrip().split()[1])
			except:
				RBPexp_absmax[i] = float(line.rstrip())
			i+=1
	for i in range(len(RBPexp_absmax)):
		rbp_exp_vec[i] = np.max([-0.99, np.min([0.99,rbp_exp_vec[i]/RBPexp_absmax[i]]) ])
		rbp_exp_vec[i+len(RBPexp_absmax)] = np.min([0.99,rbp_exp_vec[i+len(RBPexp_absmax)]/RBPexp_absmax[i]])
	return rbp_exp_vec


def make_single_table(darts_flat_fn, cis_feature_fn, rbp_fn, outfn, event_type):
	
	#### read in darts-flat output ####
	print(".. read darts-flat output")
	sp_out = read_sp_out(darts_flat_fn)
	
	#### read in rbp exp ####
	print(".. read rbp exp")
	#rbp_txt = [x for x in os.listdir(outdir) if 'RBP_tpm' in x]
	#if len(rbp_txt)>1: print("warning: more than one RBP_tpm detected in folder; only the first one is used.")
	#rbp_fn = os.path.join(outdir,rbp_txt[0])
	print(rbp_fn)
	colnames, rbp_exp_vec = read_rbp_exp(rbp_fn)
	rbp_exp_vec = RBPexp_normalizer(rbp_exp_vec, event_type)
	
	#### read in sequence features ####
	print(".. read sequence feature")
	seqFeature = read_sequence_feature(cis_feature_fn)
	colnames_2 = seqFeature.columns.values
	seqFeature_len = seqFeature.shape[1]
	colnames = np.concatenate([colnames_2, colnames])
	
	#### make the feature table ####
	print(".. make table")
	idx = [i for i in range(len(sp_out['evt'])) if sp_out['evt'][i] in seqFeature.index]
	Y = sp_out['label'][idx]
	metadata = sp_out['metadata'][idx]
	X = seqFeature.loc[[x for x in sp_out['evt'] if x in seqFeature.index]].as_matrix()	
	if X.shape[0]==0:
		raise Exception('no ID could be matched with cis- file; check your Darts-flat ID column or if you supplied the correct cis-feature file.')
	X = np.hstack([X, [rbp_exp_vec]*X.shape[0]])
	
	#### save feature table in h5 ####
	print(".. save table")
	store_fn = outfn
	with h5py.File(store_fn, 'w') as store:
		## X: n by p np-array
		## Y: n by 1 binary 0/1 np-array
		## rownames: metadata n by 8: exon_id, target, cell, experiment_id i1, s1, i2, s2
		## colnames: sequence feature and exp feature names
		store.create_dataset('X', data = X)
		store.create_dataset('Y', data = Y)
		store.create_dataset('rownames', data = metadata)
		store.create_dataset('colnames', data = np.asarray(colnames, dtype='str').astype('|S30'))


def parser( args ):
	# receive user arguments
	if args.cis is None:
		cis = config.CURRENT_CIS_PATH[args.event_type]
	else:
		cis = args.cis
	outfn = args.output
	outdir = os.path.dirname(os.path.realpath(outfn))
	darts_flat_fn = args.input
	
	# process the RBPs	
	if len(args.expr)==2:
		kallisto_1 = args.expr[0].split(',')
		kallisto_2 = args.expr[1].split(',')
		# convert kallisto transcript to gene
		logger.info('convert tx to gene TPM')
		for kal_dir in kallisto_1 + kallisto_2:
			_Darts_tx2g.parser(kal_dir)
		# extract RBP from gene
		logger.info('extract RBP')
		_Darts_rbp.parser(outdir, 
			[os.path.join(x,'gene_tpm.tsv') for x in kallisto_1],
			[os.path.join(x,'gene_tpm.tsv') for x in kallisto_2])
		rbp_fn = os.path.join(outdir, 'RBP_tpm.txt')
	else:
		rbp_fn = args.expr[0]
	
	# make the h5 store
	logger.info('make h5 feature table')
	event_type = args.event_type
	make_single_table(darts_flat_fn, cis, rbp_fn, outfn, event_type)
