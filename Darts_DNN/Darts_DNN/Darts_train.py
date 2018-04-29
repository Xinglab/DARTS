"""
Darts_DNN - Darts_train
# read in data batch and train a dnn
# Zijun Zhang
# 6.6.2017
# revised 6.9.2017: returning two random targets at a time in `read_data_batch`
# DEPRECATED as of 7.28: revised 7.25.2017: loosen filter on negative data and train a new model
# revised 7.28.2017: adapt to use Kallisto TPM
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

from keras.optimizers import RMSprop, SGD
from keras.utils import np_utils
from Darts_DNN.sklearn_nn import *

logger = logging.getLogger('Darts_DNN.train')

#### global variables ####

# number of total features
feature_num = 5922

# temporary dir
#tmp_dir = './tmp'
tmp_dir = tempfile.tempdir
logger.info('tmp_dir='+tmp_dir)

# path to all RNA-seq data
# MOVE TO USER ARGUMENTS
#target_dir_encode = '/home/zzj/scratch/Darts_deepLearning/ENCODE'
#target_dir_roadmap = '/home/zzj/scratch/Darts_deepLearning/Roadmap'

# some RBP with KD in only one cell line/didn't finish spider: 
# 58 in total
leave_out_data = ['NIP7_HepG2_ENCSR696LLZ', 'GLRX3_K562_ENCSR874DVZ', 'EIF3A_K562_ENCSR258VGD', 'GNB2L1_HepG2_ENCSR116YMU', 'HDGF_HepG2_ENCSR634KHL', 'RPLP0_HepG2_ENCSR082YGI', 'HNRNPLL_HepG2_ENCSR490DYI', 'RBM47_HepG2_ENCSR711ZJQ', 'DDX5_HepG2_ENCSR808FBR', 'TRIP6_K562_ENCSR152IWT', 'SRSF3_HepG2_ENCSR376FGR', 'KRR1_HepG2_ENCSR542ESY', 'UBE2L3_HepG2_ENCSR424JSU', 'WDR43_K562_ENCSR341PZW', 'SRSF4_K562_ENCSR697GLD', 'CCDC124_K562_ENCSR874ZLI', 'DDX51_K562_ENCSR029LGJ', 'WDR3_K562_ENCSR334BTA', 'ACO1_HepG2_ENCSR511SYK', 'ZC3H8_K562_ENCSR448JAM', 'POLR2G_K562_ENCSR936TED', 'RBM27_HepG2_ENCSR222SMI', 'AGO2_K562_ENCSR495YSS', 'PPP1R8_K562_ENCSR844QNT', 'RPS5_HepG2_ENCSR838SMC', 'RPS3_K562_ENCSR410MIQ', 'RPS2_HepG2_ENCSR667RIA', 'UTP3_HepG2_ENCSR910ECL', 'DROSHA_K562_ENCSR624XHG', 'FASTKD1_HepG2_ENCSR728BOL', 'ZNF622_HepG2_ENCSR518JXY', 'PAPOLA_K562_ENCSR368ZRP', 'ETF1_HepG2_ENCSR840QOH', 'MARK2_HepG2_ENCSR105OXX', 'HNRNPD_HepG2_ENCSR660MZN', 'MBNL1_K562_ENCSR222COT', 'SRPK2_K562_ENCSR524YXQ', 'YTHDC2_K562_ENCSR843LYF', 'SNRNP70_HepG2_ENCSR635BOO', 'KIF1C_HepG2_ENCSR182GKG', 'MAK16_K562_ENCSR517JHY', 'AGGF1_K562_ENCSR812TLY', 'RBM3_K562_ENCSR675KPR', 'RPL23A_HepG2_ENCSR706SXN', 'CPEB4_K562_ENCSR795VAK', 'TIA1_HepG2_ENCSR057GCF', 'NFX1_K562_ENCSR007XKL', 'AQR_K562_ENCSR624OUI', 'AGO3_K562_ENCSR207QGW', 'CNOT8_K562_ENCSR312HJY', 'NKRF_HepG2_ENCSR517JDK', 'CEBPZ_HepG2_ENCSR929PXS', 'PABPN1_K562_ENCSR416ZJH', 'WRN_K562_ENCSR165BCF', 'TUFM_HepG2_ENCSR459EMR', 'CALR_HepG2_ENCSR040WAK', 'NUP35_HepG2_ENCSR457WBK', 'CKAP4_HepG2_ENCSR269SJB']

# roadmap leave-out data
roadmap_leave_out = ['thymus-adipose', 'thymus-heart_right_atrium', 'thymus-adrenal_gland', 'thymus-small_intestine', 'thymus-esophagus', 'thymus-H1_cell_line', 'thymus-sigmoid_colon', 'thymus-psoas_muscle', 'thymus-H1_BMP4_derived_trophoblast_cultured_cells', 'thymus-pancreas', 'thymus-H1_BMP4_derived_mesendoderm_cultured_cells', 'thymus-heart_right_ventricle', 'thymus-liver', 'thymus-heart_aorta', 'thymus-gastric', 'thymus-bladder', 'thymus-heart_left_ventricle', 'thymus-H1_derived_mesenchymal_stem_cells', 'thymus-ovary', 'thymus-spleen', 'thymus-H1_derived_neuronal_progenitor_cultured_cells', 'thymus-lung']

# testing data record 
test_data = {
	'seen':set(),
	'X':np.empty((0,feature_num), dtype='float32'), 
	'Y':np.empty((0,2), dtype='float32'),
	'rownames':np.asarray([], dtype='str')
	}	



def construct_training_data(tar_dir, experiment_full, lite=False, annot_only=False):
	if annot_only:
		annot_evt = read_annot_evt()
	if not os.path.isfile(os.path.join(tar_dir, experiment_full, 'data.kallisto.h5')):
		logger.info('this file is not found: %s'%os.path.join(tar_dir, experiment_full, 'data.kallisto.h5'))
	with h5py.File(os.path.join(tar_dir, experiment_full, 'data.kallisto.h5'), 'r') as store:
		X = store['X'][:]
		Y = store['Y'][:]
		colnames = store['colnames'][:]
		rownames = store['rownames'][:]
	idx = [i for i in range(len(X)) if not any(np.isnan(X[i])) and
		(not annot_only or rownames[i,0] in annot_evt) and
		np.sum(np.asarray(rownames[i,4:6],dtype='int'))>20 and
		np.sum(np.asarray(rownames[i,6:8],dtype='int'))>20 and
		## we need to get rid of negative cases, so
		(Y[i]==1 or
			# min(I1, S1)>2
			(np.min(np.asarray(rownames[i,4:6],dtype='int'))>2 and
			# min(I2, S2)>2
			np.min(np.asarray(rownames[i,6:8],dtype='int'))>2)
		)
		]
	if lite:
		col_idx = []
		for i in range(len(colnames)):
			if i>2925:
				break
			if not (
					colnames[i].startswith('Yeo') or 
					colnames[i].startswith('ShortSeq') or 
					#'pssm' in colnames[i] or
					'PSSM' in colnames[i] or
					colnames[i].startswith('ISS') or  
					colnames[i].startswith('ISE') or 
					colnames[i].startswith('ESS') or 
					colnames[i].startswith('ESE')
				):
						col_idx.append(i)
	else:
		col_idx = range(X.shape[1])
	X_use = X[idx, :]
	X_use = X_use[:, col_idx]
	Y_use = Y[idx]
	return {'X':X_use, 'Y':Y_use, 'colnames':colnames[col_idx], 'rownames':rownames[idx,0]}



def write_current_pred(y_pred, fn):
	with open(fn, 'w') as fo:
		for i in range(len(test_data['rownames'])):
			fo.write('\t'.join([
				test_data['rownames'][i],
				str(test_data['Y'][i,1]),
				str(y_pred[i])
				]) + '\n')
	return


def read_data_batch_encode(target_dir_encode, n_epoch=500, verbose=True):
	experiment_list = os.listdir(target_dir_encode)
	target2exp = defaultdict(list)
	for experiment_full in experiment_list:
		if experiment_full.startswith('Control'):
			continue
		data_fn = os.path.join(target_dir_encode, experiment_full, 'data.kallisto.h5')
		if not os.path.isfile(data_fn):
			continue
		target, cell, experiment = experiment_full.split('_')
		target2exp[target].append(experiment_full)
	for epoch in range(n_epoch):
		if verbose: logger.info(">> Encode epoch "+str(epoch))
		all_targets = np.asarray([x for x in target2exp.keys() if target2exp[x][0] not in leave_out_data])
		while len(all_targets)>0:
			try:
				idx = random.sample(range(len(all_targets)),4)
			except ValueError:
				idx = range(len(all_targets))
			current_target_list = all_targets[idx]
			all_targets = np.delete(all_targets, idx)
			msg = ''
			#X_batch = np.empty((0,feature_num), dtype='float32')
			#Y_batch = np.empty((0,2), dtype='float32')
			tmp_x = []
			tmp_y = []
			for target in current_target_list:
				if len(target2exp[target])<2: continue
				for experiment_full in target2exp[target]:
					train_data = construct_training_data(target_dir_encode, experiment_full)
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
						0.1,
						12345,  ## fixing random seed so testing set is fixed
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
			Y_batch = np.array(tmp_y)
			X_batch = np.array(tmp_x)
			yield (Y_batch, X_batch, msg)
	return


def read_data_batch_roadmap(target_dir_roadmap, n_epoch=500, verbose=True):
	raw_experiment_list = os.listdir(target_dir_roadmap)
	experiment_list = []
	for exp in raw_experiment_list:
		if os.path.isfile(os.path.join(target_dir_roadmap, exp, 'data.kallisto.h5')) and\
			exp not in roadmap_leave_out:
			experiment_list.append(exp)
	for epoch in range(n_epoch):
		if verbose:  logger.info(">> Roadmap epoch "+str(epoch))
		all_targets = np.asarray([x for x in experiment_list])
		while len(all_targets)>0:
			try:
				idx = random.sample(range(len(all_targets)),1)
			except ValueError:
				idx = 0
			current_target_list = all_targets[idx] if idx!=0 else [all_targets[0]]
			all_targets = np.delete(all_targets, idx)
			msg = ''
			#X_batch = np.empty((0,feature_num), dtype='float32')
			#Y_batch = np.empty((0,2), dtype='float32')
			tmp_x = []
			tmp_y = []
			for experiment_full in current_target_list:
				train_data = construct_training_data(target_dir_roadmap, experiment_full)
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
					0.01,
					12345, ## fix random seed here
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


def train_dnn_classifier(odir='.',target_dir_encode = '/home/zzj/scratch/Darts_deepLearning/ENCODE',target_dir_roadmap = '/home/zzj/scratch/Darts_deepLearning/Roadmap'):
	# build dnn model
	clf = KerasNeuralNetClassifier(
		n_in=feature_num, 
		n_out=2,
		n_hidden=[1200, 500, 300, 200],
		dropout=[0, 0.6, 0.5, 0.3, 0.1],
		n_epoches=1,
		batch_size=1000,
		optimizer='RMSprop',
		batch_norm=True
		)
	try:
		clf.model.load_weights(os.path.join(odir,'.keras_model.h5'))
	except:
		pass
	# read each dataset iteratively
	gen_encode = read_data_batch_encode(target_dir_encode)
	if target_dir_roadmap is not None:
		gen_roadmap = read_data_batch_roadmap(target_dir_roadmap)
	index = 0
	best_val_auc = 0
	while True:
		index += 1
		try:
			Y_batch1, X_batch1, msg1 = gen_encode.next()
			if target_dir_roadmap is not None:
				Y_batch2, X_batch2, msg2 = gen_roadmap.next()
				Y_batch  = np.concatenate([Y_batch1, Y_batch2], axis=0)
				X_batch  = np.concatenate([X_batch1, X_batch2], axis=0)
				msg = msg1 + msg2
			else:
				Y_batch = Y_batch1
				X_batch = X_batch1
				msg = msg1
		except StopIteration:
			break
		
		best_auc = clf.fit(X_batch, Y_batch, 
			nb_fold=5, metric_to_use='roc', 
			random_seed=123456,
			minibatch_pos_prop=0.25
			)
		
		y_pred = clf.predict(test_data['X'])[:,1]
		auroc=metrics.roc_auc_score(test_data['Y'][:,1], y_pred)
		aupr=metrics.average_precision_score(test_data['Y'][:,1], y_pred)
		logger.info(msg)
		logger.info('auroc='+str(auroc))
		logger.info('aupr='+str(aupr))
		logger.info('val_roc='+str(best_auc))
		if auroc>best_val_auc:
			best_val_auc = auroc
			clf.model.save(os.path.join(odir,'best_val_model.h5'))
		if not index%10:
			write_current_pred(y_pred, os.path.join(tmp_dir, 'batch_{0}_pred.txt'.format(index)))
	return	

def parser( args ):
	logger.info('program starts.')
	odir = args.out_dir
	idir_list = args.train_dir
	if len(idir_list)==2:
		train_dnn_classifier(odir, idir_list[0], idir_list[1])
	else:
		train_dnn_classifier(odir, idir_list[0], None)
	logger.info('training finished.')
