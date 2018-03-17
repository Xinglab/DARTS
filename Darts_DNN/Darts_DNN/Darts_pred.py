"""
Darts_DNN - Darts_pred
# testing the dnn model on Control cell lines
# i.e. HepG2 vs. K562
# Zijun Zhang
# 6.16.2017
# revised 7.1.2017: changed to PTBP1 CRISPR
# revised 7.2.2017: changed to write prediction for a single dataset
# revised 8.8.2017: changed to kallisto TPM
# revised 2.16.2018: changed to python package
"""

import warnings  ## to filter warnings caused by Keras backward compatibility
warnings.filterwarnings('ignore')

from Darts_DNN.Darts_train import *

import logging
logger = logging.getLogger('Darts_DNN.predict')


def write_out_pred(y_pred, y_true, rownames, fn):
	with open(fn, 'w') as fout:
		fout.write("ID\tY_true\tY_pred\n")
		for i in range(len(rownames)):
			evt = rownames[i]
			this_y_pred = y_pred[i]
			this_y_true = y_true[i]
			fout.write("%s\t%f\t%f\n"%(evt, this_y_true, this_y_pred))
	return

	
def eval_continuous_label_auc(score, pred):
	idx = [i for i in range(len(score)) if score[i]>0.9 or score[i]<0.1]
	y_true = np.asarray(score[idx])
	y_true[y_true>0.5] = 1
	y_true[y_true<=0.5] = 0
	y_pred = pred[idx]
	auroc = metrics.roc_auc_score(y_true, y_pred)
	aupr = metrics.average_precision_score(y_true, y_pred)
	return auroc, aupr
	
	
def construct_data(fn, covg_filter=20):
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

def parser(args):
	# receive the arguments
	fn = args.input
	if args.model is None:
		model_fn = os.path.join(os.path.expanduser('~'), '.darts', 'trainedParam', 'encode_plus_roadmap', 'best_testing_model.h5')
	else:
		model_fn = args.model
	out_fn = args.output
	logger.info("file="+fn)
	logger.info("model="+model_fn)
	# build model and load weights
	clf = KerasNeuralNetClassifier(
		n_in=5922, 
		n_out=2,
		n_hidden=[1200, 500, 300, 200],
		dropout=[0, 0.6, 0.5, 0.3, 0.1],
		n_epoches=1,
		batch_size=400,
		optimizer='RMSprop'
		)
	
	clf.model.load_weights(model_fn)
	
	# make predictions
	data = construct_data(fn, covg_filter=0)
	pos = np.sum(data['Y']>0.9)
	neg = np.sum(data['Y']<0.1)
	y_pred = clf.predict(data['X'])[:,1]
	auroc, aupr = eval_continuous_label_auc(data['Y'], y_pred)
	
	# write predictions and log the accuracy
	write_out_pred(y_pred, data['Y'], data['rownames'], out_fn)
	logger.info("pos="+str(pos))
	logger.info("neg="+str(neg))
	logger.info("AUROC="+str(auroc))
	logger.info("AUPR="+str(aupr))
