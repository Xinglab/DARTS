# -*- coding: UTF-8 -*-

"""
Darts_DNN - predict

Predict probability of differential splicing between two conditions
"""

import sys
import warnings  ## to filter warnings caused by Keras backward compatibility
warnings.filterwarnings('ignore')

#from .Darts_train import *
from .sklearn_nn import *
from . import Darts_tx2g as _Darts_tx2g
from . import Darts_rbpTable as _Darts_rbp
from .utils import *
from . import config

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
	

def parser(args):
	# receive the arguments
	input_fn = args.input
	if args.model is None:
		model_fn = config.DEFAULT_MODEL_PATH[args.event_type]
	else:
		model_fn = args.model
	out_fn = args.output
	outdir = os.path.dirname(os.path.realpath(out_fn))
	logger.info("empirical_data="+input_fn)
	logger.info("model="+model_fn)
	# build model and load weights
	clf = KerasNeuralNetClassifier(
		**config.CURRENT_ARCHITECTURE[args.event_type]
		)
	
	clf.model.load_weights(model_fn)
	
	# make predictions
	if input_fn.endswith('h5'):
		data = construct_data_from_h5(input_fn, covg_filter=0)
	else:
		if not args.expr:
			logger.info('when input is a .txt label file, must set the expression file(s) by "-e" options. ')
			sys.exit(0)
		logger.info('reading cis and trans feaures for type %s'%args.event_type)
		if len(args.expr) == 2:
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
			geneExp_fn = os.path.join(outdir, 'RBP_tpm.txt')
		elif len(args.expr) == 1:
			geneExp_fn = args.expr[0]
		else:
			raise Exception('Two expression folders should be provided if input file is a .txt.')
		geneExp_absmax = read_geneExp_absmax(config.CURRENT_TRANS_PATH[args.event_type])
		logger.info('read cis-features')
		seqFeature_df = read_sequence_feature(config.CURRENT_CIS_PATH[args.event_type])
		logger.info('constructing in-memory feature matrix')
		data = construct_training_data_from_label(input_fn, geneExp_fn, seqFeature_df, geneExp_absmax)
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
