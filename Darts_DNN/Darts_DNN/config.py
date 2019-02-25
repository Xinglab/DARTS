# -*- coding: UTF-8 -*-

"""Here are general configurations for the Darts_DNN package, including 
version control, trained model parameter, etc.
"""

from pkg_resources import resource_filename
import os
import yaml

CURRENT_VERSION = "v0.1.0"

def substitute_home_directory(d):
	for k in d:
		if d[k].startswith('~/'):
			d[k] = os.path.join(os.path.expanduser('~'), d[k].lstrip('~/'))


## Predict and Build-feature configurations
TRANS_TOTAL_NUM = 1498
RBP_GENE_LIST_PATH = resource_filename('Darts_DNN.resources', 'rbp_gene_list.txt')
T2G_FILE_PATH = resource_filename('Darts_DNN.resources', 'human_t2g.txt')

## Train configurations
MAX_TRAIN_PATIENCE = 100

## Get_data configuations
DOWNLOAD_CONFIG = yaml.safe_load(open( resource_filename('Darts_DNN.resources', 'download.yaml'), 'r'))
CURRENT_AVAILABLE_DATA = DOWNLOAD_CONFIG[CURRENT_VERSION]

## DNN Model configurations
ARCH_CONFIG = yaml.safe_load(open( resource_filename('Darts_DNN.resources', 'architecture.yaml') ,'r'))
# for future versions, enable 'softlinks' by just providing a string
CURRENT_ARCHITECTURE = ARCH_CONFIG[ ARCH_CONFIG[CURRENT_VERSION] ] if isinstance(ARCH_CONFIG[CURRENT_VERSION], str) else ARCH_CONFIG[CURRENT_VERSION]


CIS_CONFIG = yaml.safe_load(open( resource_filename('Darts_DNN.resources', 'cis_features.yaml') ,'r'))
CURRENT_CIS_PATH = CIS_CONFIG[ CIS_CONFIG[CURRENT_VERSION] ] if isinstance(CIS_CONFIG[CURRENT_VERSION], str) else CIS_CONFIG[CURRENT_VERSION]
substitute_home_directory(CURRENT_CIS_PATH)

TRANS_CONFIG = yaml.safe_load(open( resource_filename('Darts_DNN.resources', 'trans_features.yaml') ,'r'))
CURRENT_TRANS_PATH =  TRANS_CONFIG[ TRANS_CONFIG[CURRENT_VERSION] ] if isinstance(TRANS_CONFIG[CURRENT_VERSION], str) else TRANS_CONFIG[CURRENT_VERSION]
substitute_home_directory(CURRENT_TRANS_PATH)

PARAM_CONFIG = yaml.safe_load(open( resource_filename('Darts_DNN.resources', 'trained_model_param.yaml') ,'r'))
DEFAULT_MODEL_PATH = PARAM_CONFIG[CURRENT_VERSION]
substitute_home_directory(DEFAULT_MODEL_PATH)
