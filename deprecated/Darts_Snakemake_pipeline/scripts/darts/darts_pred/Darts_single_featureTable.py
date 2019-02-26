# make the sequence and RBP features for control data
# Zijun Zhang
# 6.16.2017
# revised 7.1.2017: add a command-line argument parser
# revised 8.7.2017: change to use kallisto est.


import os
import sys
import pandas as pd
import numpy as np
import h5py
from collections import defaultdict



def read_sequence_feature():
	data=pd.read_hdf('/u/nobackup/yxing/NOBACKUP/frankwoe/ENCODE/sequence_features/ENCODE_sequenceFeature_absmax_normalized.h5')
	return data


def read_sp_out(sp_fn, sig_only=False):
	if not os.path.isfile(sp_fn):
		raise Exception("spider result not found")
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
	metadata = np.array([[x, 'Control', 'HepG2-K562', '', data.I1[x], data.S1[x], data.I2[x], data.S2[x]] for x in evt], dtype=str).astype('|S60')
	return {'evt':evt, 'label':label, 'metadata':metadata}


def read_rbp_exp(rbp_fn):
	colnames = np.empty((1498*2), dtype='str').astype('|S30')
	rbp_exp_vec = np.empty((1498*2))
	with open(rbp_fn, 'r') as fin:
		idx = 0
		firstline = True
		for line in fin:
			if firstline:
				firstline=False
				continue
			ele = line.rstrip().split()
			colnames[idx] = 'KD_' + ele[0]
			colnames[idx+1498] = 'Con_' + ele[0]
			rbp_exp_vec[idx] = float(ele[1])
			rbp_exp_vec[idx+1498] = float(ele[2])
			idx += 1
	return colnames, rbp_exp_vec

def sequenceFeature_normalizer(seqFeature):
	sequencefeature_absmax = None
	absmax = np.empty((seqFeature.shape[1],))
	if sequencefeature_absmax is None:
		if not os.path.isfile('sequencefeature_absmax.txt'):
			for i in range(seqFeature.shape[1]):
				maxval = np.max(np.abs(seqFeature.iloc[:,i]))
				if maxval==0:
					print("no variance in seqFeature column: "+str(i))
					maxval=1
				absmax[i] = maxval
				#if len(set(seqFeature.iloc[:,i]))==2: # if is binary feature
				#	absmax[i] = 1.0
				#else:
				#	absmax[i] = np.percentile(np.abs(seqFeature.iloc[:,i]), 99)
			print(".. get all absmax")
			with open('sequencefeature_absmax.txt', 'w') as fout:
				fout.write('\n'.join([str(x) for x in absmax]))
			
		else:
			with open('sequencefeature_absmax.txt', 'r') as fin:
				i=0
				for line in fin:
					absmax[i] = float(line.rstrip())
					i+=1
		sequencefeature_absmax = absmax
	for i in range(seqFeature.shape[0]):
		if not i%1000:
			print(i)
		seqFeature.iloc[i,:] = seqFeature.iloc[i,:]/sequencefeature_absmax
	
	return seqFeature


def RBPexp_normalizer(rbp_exp_vec):
	#global RBPexp_absmax

	src_dir = os.path.dirname(os.path.realpath(__file__))
	rbp_max_fn = os.path.join(src_dir, 'RBPexp_absmax.kallisto.bak20170916.txt')
	#rbp_max_fn = os.path.join(src_dir, 'RBPexp_absmax.kallisto.bak.encode_norm.txt')
	RBPexp_absmax = np.zeros((1498))
	with open(rbp_max_fn, 'r') as fin:
		i=0
		for line in fin:
			#absmax[i] = float(line.rstrip())
			try:
				RBPexp_absmax[i] = float(line.rstrip().split()[1])
			except:
				RBPexp_absmax[i] = float(line.rstrip())
			i+=1
	for i in range(len(RBPexp_absmax)):
		rbp_exp_vec[i] = np.max([-0.99, np.min([0.99,rbp_exp_vec[i]/RBPexp_absmax[i]]) ])
		rbp_exp_vec[i+len(RBPexp_absmax)] = np.min([0.99,rbp_exp_vec[i+len(RBPexp_absmax)]/RBPexp_absmax[i]])
	return rbp_exp_vec


def make_single_table(sp_dir, rbp_dir, outdir):
	
	#### read in spider output ####
	print(".. read spider output")
	sp_fn = os.path.join(sp_dir,'Sp_out.txt')
	sp_out = read_sp_out(sp_fn)
	
	#### read in rbp exp ####
	print(".. read rbp exp")
	rbp_txt = [x for x in os.listdir(rbp_dir) if 'RBP_tpm' in x]
	if len(rbp_txt)>1: print("warning: more than one RBP_tpm detected in folder.")
	rbp_fn = os.path.join(rbp_dir,rbp_txt[0])
	print(rbp_fn)
	colnames, rbp_exp_vec = read_rbp_exp(rbp_fn)
	rbp_exp_vec = RBPexp_normalizer(rbp_exp_vec)
	
	#### read in sequence features ####
	print(".. read sequence feature")
	seqFeature = read_sequence_feature()
	#print ".. normalizing"
	#seqFeature = sequenceFeature_normalizer(seqFeature)
	colnames_2 = seqFeature.columns.values
	seqFeature_len = seqFeature.shape[1]
	colnames = np.concatenate([colnames_2, colnames])
	
	#### make the feature table ####
	print(".. make table")
	idx = [i for i in range(len(sp_out['evt'])) if sp_out['evt'][i] in seqFeature.index]
	Y = sp_out['label'][idx]
	metadata = sp_out['metadata'][idx]
	X = seqFeature.loc[[x for x in sp_out['evt'] if x in seqFeature.index]].as_matrix()	
	X = np.hstack([X, [rbp_exp_vec]*X.shape[0]])
	
	#### save feature table in h5 ####
	print(".. save table")
	store_fn = os.path.join(outdir, 'data.h5')
	with h5py.File(store_fn, 'w') as store:
		## X: n by p np-array
		## Y: n by 1 binary 0/1 np-array
		## rownames: metadata n by 8: exon_id, target, cell, experiment_id i1, s1, i2, s2
		## colnames: sequence feature and exp feature names
		store.create_dataset('X', data = X)
		store.create_dataset('Y', data = Y)
		store.create_dataset('rownames', data = metadata)
		store.create_dataset('colnames', data = np.asarray(colnames, dtype='str').astype('|S30'))

def parser(sp_dir, rbp_dir, outdir):
	make_single_table(sp_dir, rbp_dir, outdir)

if __name__ == '__main__':
	#make_single_table('/home/zzj/scratch/ENCODE/Control_HepG2_v_K562/Full_merge')
	#make_single_table('/u/home/f/frankwoe/nobackup/ENCODE/Control_HepG2_v_K562/Full_merge')
	#make_single_table('/u/home/f/frankwoe/nobackup/ENCODE/Control_HepG2_v_K562/Down25_merge')
	sp_dir, rbp_dir, outdir = sys.argv[1:4]
	make_single_table(sp_dir, rbp_dir, outdir)
