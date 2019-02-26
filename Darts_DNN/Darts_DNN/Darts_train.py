# -*- coding: UTF-8 -*-

"""
Darts_DNN - Darts_train

read in data batches and train a dnn
"""

import warnings  ## to filter warnings caused by Keras backward compatibility
warnings.filterwarnings('ignore')

import sys
import os
import h5py
from collections import defaultdict
import logging
import datetime
import tempfile
import random
from keras.utils import np_utils
from .sklearn_nn import *
from . import config

logger = logging.getLogger('Darts_DNN.train')
# temporary dir
tmp_dir = tempfile.tempdir
logger.info('tmp_dir='+tmp_dir)


def write_current_pred(test_data, y_pred, fn):
	with open(fn, 'w') as fo:
		for i in range(len(test_data['rownames'])):
			fo.write('\t'.join([
				test_data['rownames'][i],
				str(test_data['Y'][i,1]),
				str(y_pred[i])
				]) + '\n')
	return


def read_data_batch(train_fp, test_data, data_source_id, seqFeature_df, geneExp_absmax, n_epoch=500, verbose=True):
	''' given a train_file, return a generator
	'''
	train_list = read_train_list(train_fp)
	for epoch in range(n_epoch):
		if verbose:  logger.info(">> {} epoch {}".format(data_source_id, epoch))
		all_targets = np.asarray([x for x in train_list.index])
		while len(all_targets)>0:
			try:
				idx = random.sample(range(len(all_targets)), 4)
			except ValueError:  # not enough for sampling
				idx = np.arange(len(all_targets))
			current_target_list = all_targets[idx]
			all_targets = np.delete(all_targets, idx)
			msg = ''
			tmp_x = []
			tmp_y = []
			for experiment_full in current_target_list:
				ID = experiment_full
				if 'label_fn' in train_list:   
					## train_list is txt file, and construct feature matrix in-memory
					label_fn = train_list.loc[ID].label_fn
					geneExp_fn = train_list.loc[ID].geneExp_fn
					train_data = construct_training_data_from_label(
						label_fn, 
						geneExp_fn, 
						seqFeature_df, 
						geneExp_absmax, 
						ID=ID,
						in_training_phase=True)
				elif 'input_fn' in train_list:
					## train_list is  h5 file, and read in pre-built feature matrix from disk
					input_fn = train_list.loc[ID].input_fn
					train_data = construct_training_data_from_h5(
						input_fn,
						in_training_phase=True)
				else:
					raise Exception('input_filelist format not understood. Check online documentation for examples.')
				if verbose: 
					print("... "+experiment_full)
					print("... pos={0}, neg={1}".format(
						np.sum(train_data['Y']==1),
						np.sum(train_data['Y']==0),))
					msg += "\n{0}: pos={1}, neg={2}".format(
						experiment_full,
						np.sum(train_data['Y']==1),
						np.sum(train_data['Y']==0)
						)
				train_data['Y'] = np_utils.to_categorical(train_data['Y'], 2)
				X_test, Y_test, X_train, Y_train, arg_split = split_imbalance_data(
					train_data['X'], 
					train_data['Y'], 
					0.05,
					12345, ## fix random seed here for consistent testing set
					train_data['rownames']
					)
				rowname_test, rowname_train = arg_split[0]
				rowname_test = [x+'_'+experiment_full for x in rowname_test]
				if not experiment_full in test_data['seen']:
					test_data['X'] = np.concatenate([test_data['X'], X_test], axis=0)
					test_data['Y'] = np.concatenate([test_data['Y'], Y_test], axis=0)
					test_data['rownames'] = np.concatenate([test_data['rownames'],rowname_test], axis=0)
					test_data['seen'].add(experiment_full)
				#X_batch = np.concatenate([X_batch, X_train], axis=0)
				tmp_x.extend(X_train)
				#Y_batch = np.concatenate([Y_batch, Y_train], axis=0)
				tmp_y.extend(Y_train)
			X_batch = np.array(tmp_x)
			Y_batch = np.array(tmp_y)
			yield (Y_batch, X_batch, msg)
	return


def train_dnn_classifier(event_type, train_filelist, odir='.'):
	# build dnn model
	clf = KerasNeuralNetClassifier(
		**config.CURRENT_ARCHITECTURE[event_type]
		)
	try:
		clf.model.load_weights(os.path.join(odir,'.keras_model.h5'))
		logger.info("found and loaded previous intermediate parameters")
	except:
		pass
	# keep a testing data record 
	feature_num = config.CURRENT_ARCHITECTURE[event_type]['n_in']
	output_num  = config.CURRENT_ARCHITECTURE[event_type]['n_out']
	test_data = {
		'seen':set(),
		'X':np.empty((0,feature_num), dtype='float32'), 
		'Y':np.empty((0,output_num), dtype='float32'),
		'rownames':np.asarray([], dtype='str')
		}

	# read in seqFeature and geneExp_absmax for reference of different generators
	geneExp_absmax = read_geneExp_absmax(config.CURRENT_TRANS_PATH[event_type])
	seqFeature_df = read_sequence_feature(config.CURRENT_CIS_PATH[event_type])

	# prepare the generator each dataset iteratively
	gen_list = [
		read_data_batch( 
			train_filelist[i], 
			test_data, 
			"data_source_%i"%i,
			seqFeature_df,
			geneExp_absmax) 
		for i in range(len(train_filelist))
		]
	step = 0
	best_val_auc = 0
	patience = 0
	while True:
		step += 1
		tmp_X = []
		tmp_y = []
		msg = ''
		for gen in gen_list:
			try:
				Y_batch, X_batch, msg_batch = next(gen)
				tmp_X.append(X_batch)
				tmp_y.append(Y_batch)
				msg += msg_batch
			except StopIteration:
				break
		
		X = np.concatenate(tmp_X, axis=0)
		Y = np.concatenate(tmp_y, axis=0)
		best_auc = clf.fit(X, Y, 
			nb_fold=5, metric_to_use='roc', 
			random_seed=123456,
			minibatch_pos_prop=0.25
			)
		
		y_pred = clf.predict(test_data['X'])[:,1]
		auroc=metrics.roc_auc_score(test_data['Y'][:,1], y_pred)
		aupr=metrics.average_precision_score(test_data['Y'][:,1], y_pred)
		logger.info(msg)
		logger.info('auroc=%.3f'%(auroc))
		logger.info('aupr=%.3f'%(aupr))
		logger.info('val_roc=%.3f'%(best_auc))
		if auroc > best_val_auc:
			best_val_auc = auroc
			patience = 0
			clf.model.save(os.path.join(odir,'{}.trainedParam.h5'.format(event_type)))
		else:
			patience += 1
			if patience > config.MAX_TRAIN_PATIENCE:
				logger.info('early stop @ train episode %i.'%step)
				break
		if not step%10:
			write_current_pred(test_data, y_pred, os.path.join(tmp_dir, 'batch_{0}_pred.txt'.format(step)))
	return	

def parser( args ):
	# parse options
	logger.info('program starts.')
	odir = args.out_dir
	train_filelist = args.input_filelist
	event_type = args.event_type
	assert len(train_filelist)<=2, Exception('current version accepts at most 2 train file lists.')

	# training
	train_dnn_classifier(event_type, train_filelist, odir)

	logger.info('training finished.')
