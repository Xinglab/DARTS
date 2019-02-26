# testing the dnn model on leave-out-samples
# Zijun Zhang
# 6.12.2017
# revised 7.4.2017: add plotting function

from deep_neural_net_gpu import *

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn import metrics
import numpy as np

#### perturbation function ####

def permute_location(X, colnames, counter=None):
	seq_idx = []
	rbp_idx = []
	new_X = np.copy(X)
	for j in range(len(colnames)):
		if colnames[j].startswith('Con_') or colnames[j].startswith('KD_'):
			rbp_idx.append(j)
		else:
			seq_idx.append(j)
	random.shuffle(rbp_idx)
	new_X = X[:,seq_idx+rbp_idx]
	return new_X

def permute_withinRow(X, counter=None):
	idx = range(X.shape[0])
	new_X = np.copy(X)
	#for j in range(len(colnames)):
	#	if colnames[j].startswith('Con_') or colnames[j].startswith('KD_'):
	#		random.shuffle(idx)
	#		new_X[:,j] = X[:,j][idx]
	random.shuffle(idx)
	new_X[:,counter] = X[:,counter][idx]
	return new_X

def swap_location(X, counter=None):
	new_X = np.copy(X)
	new_X[:,counter[1]] = X[:,counter[0]]
	new_X[:,counter[0]] = X[:,counter[1]]
	return new_X

def in_sillico_perturbation(X, colnames, num=1):
	seq_idx = []
	rbp_idx = []
	new_X = np.copy(X)
	for j in range(len(colnames)):
		if colnames[j].startswith('Con_') or colnames[j].startswith('KD_'):
			rbp_idx.append(j)
		else:
			seq_idx.append(j)
	perturb_idx = random.sample(rbp_idx, num)
	perturb_col = colnames[perturb_idx]
	changed_value = 0.9 if random.random()>0.5 else 0.1
	new_X[:, perturb_idx] = changed_value
	return new_X, perturb_col


#perturb_num=100
#X_new, perturb_col = in_sillico_perturbation(train_data['X'], train_data['colnames'], perturb_num)
#print perturb_col
#y_pred2 = clf.predict(X_new)[:,1]
#metrics.roc_auc_score(train_data['Y'], y_pred2)

#np.mean(y_pred2[train_data['Y']==1])
#np.mean(y_pred2[train_data['Y']==0])


target_dir_roadmap75 = '/home/zzj/scratch/Darts_deepLearning/Roadmap_75'
roadmap_75_exp = os.listdir(target_dir_roadmap75)

#### build model and load weights ####
clf = KerasNeuralNetClassifier(
	n_in=5922, 
	n_out=2,
	n_hidden=[1200, 500, 300, 200],
	dropout=[0, 0.6, 0.5, 0.3, 0.1],
	n_epoches=1,
	batch_size=400,
	optimizer='RMSprop'
	)
#clf.model.load_weights('best_val_model.h5')
clf.model.load_weights('best_testing_model.h5')

tar = sys.argv[1]
tar = tar.lower()

#### testing on leave-out-data
print '\t'.join(['ID', '#pos', '#neg', 'AUROC', 'AUPR'])
counter = 0
#gen = leave_out_data if tar == 'encode' else roadmap_leave_out
if tar == 'encode':
	gen = leave_out_data
elif tar == 'roadmap':
	gen = roadmap_leave_out
elif tar == 'roadmap75':
	gen = roadmap_75_exp
else:
	raise Exception("i don't understand set %s"%tar)

for experiment_full in gen:
#for experiment_full in os.listdir(target_dir):
	if experiment_full != 'thymus-liver':
		continue
	counter += 1
	try:
		if tar.lower() == 'encode':
			train_data = construct_training_data(target_dir_encode, experiment_full)
		elif tar.lower() == 'roadmap':
			train_data = construct_training_data(target_dir_roadmap, experiment_full)
		elif tar.lower() == 'roadmap75':
			train_data = construct_training_data(target_dir_roadmap75, experiment_full)
		else:
			raise Exception("i don't understand set '%s'"%tar)
	except:
		continue
	pos = np.sum(train_data['Y']==1)
	neg = np.sum(train_data['Y']==0)
	y_pred = clf.predict(train_data['X'])[:,1]
	auroc = metrics.roc_auc_score(train_data['Y'], y_pred)
	aupr = metrics.average_precision_score(train_data['Y'], y_pred)
	print >> sys.stdout, '\t'.join([experiment_full, str(pos), str(neg), str(auroc), str(aupr)])
	
	#plotting
	#fpr, tpr, thresh = metrics.roc_curve(train_data['Y'], y_pred)
	#auc = metrics.roc_auc_score(train_data['Y'], y_pred)
	#plt.plot(fpr,tpr,label=target+' '+cell+", auc="+str(int(auc*100)/100.))

#plt.legend(loc=0)
#plt.savefig('leave_out_auroc.pdf')
