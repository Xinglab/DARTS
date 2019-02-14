# -*- coding: UTF-8 -*-

'''
Darts_DNN - utils
training and prediction utilities

Below is original header:
# utilities for machine learning stuff
# Zijun Zhang, 1.19.2017
# last revised: 2.22.2017
# revised 6.6.2017: reset random seed to None after split fold
# revised 6.6.2017: fixed a bug in `split_data_fold` when auxillary arguments are provided
# revised 6.7.2017: provide a proportion for minibatch construction in `construct_balanced_minibatch`
# revised 6.7.2017: moved `split_imbalance_data` into utils
'''

import os
import random
import pandas as pd
import numpy as np
import keras.backend as K
from sklearn import metrics
import h5py
import scipy.stats
from . import config

def split_data_fold(X, Y, n_fold=5, fold_id=0, shuffle=True, random_seed=123, *args):
	""" split data X, Y into n_fold and return one fold versue the rest
		X: numpy ndarray,
		Y: list/numpy array
		n_fold: total number of fold to split
		fold_id: fold to be returned
		shuffle: randomly re-order if turned on
		*args: other axillary attributes to be split
	"""
	assert fold_id < n_fold
	random.seed(random_seed)
	idx = range(X.shape[0])
	if shuffle: 
		for n in range(100): random.shuffle(idx)
	_X = np.asarray(X[idx,], dtype=K.floatx())
	_Y = np.asarray([Y[i] for i in idx], dtype='float32')
	fold_size = len(idx) // n_fold
	target_set = set(range(fold_id*fold_size, (fold_id+1)*fold_size))
	remain_set = set(range(len(idx))) - target_set
	X_target = np.asarray(_X[list(target_set), ], dtype=K.floatx())
	X_remain = np.asarray(_X[list(remain_set), ], dtype=K.floatx())
	Y_target = np.asarray(_Y[list(target_set)], dtype='float32')
	Y_remain = np.asarray(_Y[list(remain_set)], dtype='float32')
	arg_split = []
	for arg in args:
		arg = np.asarray(arg)[idx]
		arg_target = arg[list(target_set)]
		arg_remain = arg[list(remain_set)]
		arg_split.append([arg_target, arg_remain])
	
	random.seed(None)
	if len(arg_split)>0:	
		return X_target, Y_target, X_remain, Y_remain, arg_split
	else:
		return X_target, Y_target, X_remain, Y_remain


def split_imbalance_data(X, Y, prop_test=0.1, random_seed=123456, *args):	
	""" randomly choose the given proportion of positive labels out for testing, with 
	    matching amount of negative labels
		X: np-array, features
		Y: np-array
		prop_test: float, proportion of positive labels for testing
		random_seed: int
		*args: other axillary attributes to be split
	"""
	random.seed(random_seed)
	r0_index = split_y_by_label(Y, 0)
	r1_index = split_y_by_label(Y, 1)
	idx = range(X.shape[0])
	
	for n in range(10): random.shuffle(r0_index)
	for n in range(10): random.shuffle(r1_index)
	
	_X = np.asarray(X, dtype='float32')
	_Y = np.asarray(Y, dtype='float32')
	
	fold_size = int(len(r1_index)*prop_test)
	target_set = set(np.concatenate([r0_index[0:fold_size], r1_index[0:fold_size]], axis=0))
	remain_set = set(range(len(idx))) - target_set
	X_target = np.asarray(_X[list(target_set), ], dtype='float32')
	X_remain = np.asarray(_X[list(remain_set), ], dtype='float32')
	Y_target = np.asarray(_Y[list(target_set)], dtype='float32')
	Y_remain = np.asarray(_Y[list(remain_set)], dtype='float32')
	arg_split = []
	for arg in args:
		arg = np.asarray(arg)
		arg_target = arg[list(target_set)]
		arg_remain = arg[list(remain_set)]
		arg_split.append([arg_target, arg_remain])
	
	random.seed(None)
	if len(arg_split)>0:	
		return X_target, Y_target, X_remain, Y_remain, arg_split
	else:
		return X_target, Y_target, X_remain, Y_remain
	return
		
	
def split_data_tier(X, Y, R, T=2, percents=None, shuffle=True, verbose=True):
	""" based on the median/percentile of R (auxillary factor),
		split X and Y seperately for Y=0 and Y=1
		X: numpy array
		Y: list/numpy array
		R: list/numpy array, auxillary factor
	"""
	if percents is None: percents = np.concatenate([np.arange(0., 1., step=1./T)[1:],[1.]]) * 100.
	R0_percentiles = np.percentile([R[i] for i in range(len(Y)) if Y[i]==0], percents)
	R1_percentiles = np.percentile([R[i] for i in range(len(Y)) if Y[i]==1], percents)
	if verbose: print [R0_percentiles, R1_percentiles]
	R0_index = []
	R1_index = []
	prev_percent = -np.inf
	for percent in R0_percentiles:
		R0_index.append([i for i in range(len(R)) if Y[i]==0 and R[i]>prev_percent and R[i]<=percent])
		prev_percent = percent
	prev_percent = -np.inf
	for percent in R1_percentiles:
		R1_index.append([i for i in range(len(R)) if Y[i]==1 and R[i]>prev_percent and R[i]<=percent])
		prev_percent = percent
	if shuffle:
		for j in range(len(R0_index)): np.random.shuffle(R0_index[j])
		for j in range(len(R1_index)): np.random.shuffle(R1_index[j])
	return R0_index, R1_index


def construct_balanced_minibatch(index, R0_idx, R1_idx, batch_size, pos_prop=0.5, shuffle=True):
	""" construct a mini-batch from R[01]_index to be balanced (50% 0, 50% 1)
		index: int, index of mini-batch
		R[01]_idx: a list/numpy array for the two classes
		batch_size: int
		pos_prot: float, the proportion of positive labels in a batch
	"""
	b1_size = int(batch_size*pos_prop)
	b0_size = batch_size - b1_size
	n0 = len(R0_idx) 
	n1 = len(R1_idx)
	if (index*b0_size%n0) < ((index+1)*b0_size%n0):
		R0_sam = R0_idx[(index*b0_size%n0):((index+1)*b0_size%n0)]
	else:
		R0_sam = R0_idx[(index*b0_size%n0):] + R0_idx[:((index+1)*b0_size%n0)]
	if (index*b1_size%n1) < ((index+1)*b1_size%n1):
		R1_sam = R1_idx[(index*b1_size%n1):((index+1)*b1_size%n1)]
	else:
		R1_sam = R1_idx[(index*b1_size%n1):] + R1_idx[:((index+1)*b1_size%n1)]
	sam = R0_sam + R1_sam
	if shuffle: random.shuffle(sam)
	return sam


def split_y_by_label(y, label):
	""" a utility function to automatically 
	find the indices of y for given a label
	"""
	unique_labels = np.unique(y)
	if np.min(unique_labels)==-1: label = label*2-1
	if len(y.shape)>1:
		return [i for i in range(len(y)) if y[i,1]==label]
	else:
		return [i for i in range(len(y)) if y[i]==label]


def train_keras_balanced_model(model, x_train, y_train, x_val, y_val, 
		sample_weight=None, n_epoch=100, batch_size=64, minibatch_pos_prop=0.5, 
		save_model_name='my_model.h5', metric_to_use='roc'):
	""" train keras model by `train_on_batch` while using the user-defined
	biased-proportion of minibatches and custom early-stop surveillance
	"""
	r0_index = split_y_by_label(y_train, 0)
	r1_index = split_y_by_label(y_train, 1)
	
	n_R1_batches = len(r1_index)//(batch_size*minibatch_pos_prop)
	n_R0_batches = len(r0_index)//(batch_size*(1-minibatch_pos_prop))
	n_train_batches = int(max(n_R1_batches, n_R0_batches))
	epoch = 0
	index = 0
	best_iter = 0
	best_auc = -np.inf if metric_to_use!='loss' else np.inf
	#metric_to_use = 'pr'
	
	patience = 1000
	patience_increase = 2
	improvement_threshold = 1.01 if metric_to_use!='loss' else 0.95
	#validation_frequency = min(n_R1_batches, patience //2)
	validation_frequency = 10
	
	#print x_train.shape
	done_looping = False
	while (epoch < n_epoch) and (not done_looping):
		epoch += 1
		for index in range(n_train_batches):
			iter = (epoch - 1) * n_train_batches + index
			sam_idx = construct_balanced_minibatch(index, r0_index, r1_index, batch_size, pos_prop=minibatch_pos_prop)
			if sample_weight is not None:
				model.train_on_batch(
					x_train[sam_idx,:], 
					y_train[sam_idx],
					sample_weight=sample_weight[sam_idx])
			else:
				model.train_on_batch(
					x_train[sam_idx,:], 
					y_train[sam_idx])
			#if iter % 10 == 0:
			#	print('training @ iter = ' + str(iter))
			print('training @ iter = ' + str(iter))
			if (not iter % validation_frequency) or index==n_train_batches-1: 
				this_y_val_hat = model.predict(x_val)
				this_val_eval = model.evaluate(x_val, y_val)
				this_train_eval = model.evaluate(x_train, y_train)
				if len(y_val.shape)>1:
					this_roc = metrics.roc_auc_score(y_val[:,1], this_y_val_hat[:,1])
					this_pr = metrics.average_precision_score(y_val[:,1], this_y_val_hat[:,1])
				else:
					this_roc = metrics.roc_auc_score(y_val, this_y_val_hat)
					this_pr = metrics.average_precision_score(y_val, this_y_val_hat)
				this_metrics = {'acc':this_val_eval[1], 'loss':this_val_eval[0],
					'roc':this_roc, 'pr':this_pr}
				this_auc = this_metrics[metric_to_use]
				print('epoch %i, minibatch %i/%i, train loss %.3f, validation loss %.3f, auc %.3f %%' %
					(epoch, index + 1, n_train_batches,
					this_train_eval[0],
					this_val_eval[0], this_auc))
				if (this_auc > best_auc and metric_to_use!='loss') or (this_auc < best_auc and metric_to_use=='loss'):
					if metric_to_use != 'loss':
						if this_auc > best_auc * improvement_threshold:
							patience = max(patience, iter*patience_increase)
					else:
						if this_auc < best_auc * improvement_threshold:
							patience = max(patience, iter*patience_increase)
					best_auc = this_auc
					best_iter = iter
					model.save(save_model_name)
					print "-"*50
					print epoch
					print ("     epoch %i, minibatch %i/%i") % (epoch, index+1, n_train_batches)
					print "best auc = "+str(best_auc)
					print "-"*50
			if patience <= iter:
				done_looping=True
				break
		#epoch += 1
	model.load_weights(save_model_name)
	print "Optimization complete. Best auc="+str(best_auc)+'@ iter ' +str(best_iter)
	return best_auc


def construct_data_from_h5(fn, covg_filter=20):
	with h5py.File(fn, 'r') as store:
		X = store['X'][:]
		Y = store['Y'][:]
		colnames = store['colnames'][:]
		rownames = store['rownames'][:]
	idx = [i for i in range(len(X)) if not any(np.isnan(X[i])) and
		np.sum(np.asarray([int(float(x)) for a in rownames[i,4:6] for x in a.split(',') ]))>covg_filter and
		np.sum(np.asarray([int(float(x)) for a in rownames[i,6:8] for a in a.split(',') ]))>covg_filter #and
		## we need to get rid of negative cases, so
		## THIS WAS USED IN TRAINING; BUT DEPRECATED IN PREDICTION
		#(Y[i]>0.9 or
		#	# min(I1, S1)>2
		#	(np.min(np.asarray(rownames[i,4:6],dtype='int'))>2 and
		#	# min(I2, S2)>2
		#	np.min(np.asarray(rownames[i,6:8],dtype='int'))>2)
		#)
		]
	col_idx = range(X.shape[1])
	X_use = X[idx, :]
	X_use = X_use[:, col_idx]
	Y_use = Y[idx]
	return {'X':X_use, 'Y':Y_use, 'rownames':rownames[idx,0]}


def construct_training_data_from_label(label_fn, geneExp_fn, seqFeature_df, geneExp_absmax):
	if not os.path.isfile(label_fn):
		raise Exception('this file is not found: %s'%label_fn)
	if not os.path.isfile(geneExp_fn):
		raise Exception('this file is not found: %s'%geneExp_fn)
	# read in two sets of labelled events
	# from the DARTS BHT-flat output
	pos, neg = read_label_fn(label_fn)
	pos_neg = pos + neg
	# read in the gene expression for this
	# comparison, and make it a matrix of
	# len(pos)+len(neg) rows
	geneExp, geneExp_colnames = read_geneExp(geneExp_fn, geneExp_absmax, nrow=len(pos)+len(neg) )
	# extract  seqfeatures
	seqFeature = seqFeature_df.loc[pos+neg]
	assert seqFeature.shape[0] == len(pos) + len(neg)	
	# concat
	pos_neg_idx = [i for i in range(len(pos_neg)) if pos_neg[i] in seqFeature_df.index ]
	X_use = np.concatenate([seqFeature, geneExp], axis=1)[pos_neg_idx]
	Y_use = np.asarray([1.]*len(pos) + [0.]*len(neg))[pos_neg_idx]
	rownames = np.asarray([ x for x in pos + neg ])[pos_neg_idx]
	colnames = np.concatenate([seqFeature_df.columns.values, geneExp_colnames], axis=0)
	return {'X':X_use, 'Y':Y_use, 'colnames':colnames, 'rownames':rownames}


def read_sequence_feature(fn=None):
	if not os.path.isfile(fn):
		raise Exception('cis feature file not found: %s'%fn)
	if fn.endswith('h5') or fn.endswith('hdf5'):
		data=pd.read_hdf(fn)
	elif fn.endswith('txt'):
		data=pd.read_table(fn)
		data.index = data.ID
		data.drop(columns='ID', inplace=True)
	return data


def read_label_fn(label_fn, significance=0.1, min_read_cov=20):
	pos = set()
	neg = set()
	with open(label_fn, 'r') as f:
		firstline = True
		for line in f:
			ele = line.strip().split()
			if firstline:
				header = {ele[i]:i for i in range(len(ele))}
				firstline = False
				continue
			if ele[header['post_pr']] == "NA":
				continue
			post_pr = float(ele[header['post_pr']])
			eid = ele[header['ID']]
			I1 = sum([int(x) for x in ele[header['I1']].split(',') ])
			I2 = sum([int(x) for x in ele[header['I2']].split(',') ])
			S1 = sum([int(x) for x in ele[header['S1']].split(',') ])
			S2 = sum([int(x) for x in ele[header['S2']].split(',') ])
			if I1+S1 <= min_read_cov or I2+S2 <= min_read_cov:
				continue 
			if post_pr > 1-significance:
				pos.add(eid)
			elif post_pr < significance and min(I1,S1)>2 and min(I2,S2)>2:
				neg.add(eid)
	return list(pos), list(neg)


def read_label_fn_ad_hoc(label_fn, significance=0.1, min_read_cov=20):
	pos = set()
	neg = set()
	with open(label_fn, 'r') as f:
		firstline = True
		for line in f:
			ele = line.strip().split()
			if firstline:
				header = {ele[i]:i for i in range(len(ele))}
				firstline = False
				continue
			post_pr = float(ele[header['post_pr']])
			eid = ele[header['ID']]
			I1 = float(ele[header['I1']])
			I2 = float(ele[header['I2']])
			S1 = float(ele[header['S1']])
			S2 = float(ele[header['S2']])
			inc_len = float(ele[header['inc_len']])
			skp_len = float(ele[header['skp_len']])
			psi1 = I1/inc_len / (I1/inc_len + S1/skp_len)
			psi2 = I2/inc_len / (I2/inc_len + S2/skp_len)
			if I1+S1 <= min_read_cov or I2+S2 <= min_read_cov:
				continue 
			if (I1+S1) <= 20 and (I2+S2) <= 20:
				continue

			if abs(psi1 - psi2) > 0.15:
				pos.add(eid)
			elif abs(psi1 - psi2) < 0.015 and min(I1,S1)>2 and min(I2,S2)>2:
				neg.add(eid)
	return list(pos), list(neg)


def read_train_list(fn):
	df = pd.read_table(fn, index_col = 0)
	return df


def read_geneExp_absmax(fn):
	df = pd.read_table(fn, index_col=0, header=None).values.reshape((config.TRANS_TOTAL_NUM))
	return df


def read_geneExp(geneExp_fn, geneExp_absmax, nrow=1):
	df = pd.read_table(geneExp_fn, index_col=0).transpose()
	colnames = [prefix+'-'+suffix for prefix in df.index for suffix in df.columns.values]
	vec = df.values.flatten() / np.tile(geneExp_absmax,2 )
	n_element = len(vec)
	# mat = np.repeat(vec, nrow).reshape((n_element, nrow)).transpose()
	mat = np.tile(vec, (nrow, 1))
	return mat, colnames


