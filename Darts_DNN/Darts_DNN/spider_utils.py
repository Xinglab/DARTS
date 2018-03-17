# utilities for machine learning stuff
# Zijun Zhang, 1.19.2017
# last revised: 2.22.2017
# revised 6.6.2017: reset random seed to None after split fold
# revised 6.6.2017: fixed a bug in `split_data_fold` when auxillary arguments are provided
# revised 6.7.2017: provide a proportion for minibatch construction in `construct_balanced_minibatch`
# revised 6.7.2017: moved `split_imbalance_data` into utils

import random
import numpy as np
import theano
from sklearn import metrics
import scipy.stats

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
	_X = np.asarray(X[idx,], dtype=theano.config.floatX)
	_Y = np.asarray([Y[i] for i in idx], dtype='float32')
	fold_size = len(idx) // n_fold
	target_set = set(range(fold_id*fold_size, (fold_id+1)*fold_size))
	remain_set = set(range(len(idx))) - target_set
	X_target = np.asarray(_X[list(target_set), ], dtype=theano.config.floatX)
	X_remain = np.asarray(_X[list(remain_set), ], dtype=theano.config.floatX)
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


def co_train_keras_model(model1, x1_train, y1_train, x1_val, y1_val, batch_size1,
						model2, x2_train, y2_train, x2_val, y2_val, batch_size2,
						n_epoch=20, 
						save_model1_name='my_model_1.h5', save_model2_name='my_model_2.h5',
						w2_train = None,
						ensemble=True,
						fix_shared_layer=True,
						enable_sharing=False,
						shuffle=True
						):
	"""
	model1: baseline model
	model2: axillary model
	fix_shared_layer: all shared layer weights fixed at model1
	enable_sharing: simultaneously update shared layer in two models
	"""
	def compute_train_batches(y, batch_size):
		r0_index = [i for i in range(len(y)) if y[i]==-1]
		r1_index = [i for i in range(len(y)) if y[i]==1]
		n_R1_batches = len(r1_index)//(batch_size/2)
		n_R0_batches = len(r0_index)//(batch_size/2)
		n_train_batches = max(n_R1_batches, n_R0_batches)
		return n_train_batches, r0_index, r1_index
	
	n_train1_batches, r10_index, r11_index = compute_train_batches(y1_train, batch_size1)
	n_train2_batches, r20_index, r21_index = compute_train_batches(y2_train, batch_size2)
	n_train_batches = max(n_train1_batches, n_train2_batches)
	
	layer_dict1 = dict([(layer.name, layer) for layer in model1.layers])
	layer_dict2 = dict([(layer.name, layer) for layer in model2.layers])
	#fx_W, fx_b = layer_dict1['fx'].get_weights()
	shared_layer = {}
	for layer_name in layer_dict2:
		if layer_name in layer_dict1:
			shared_layer[layer_name] = layer_dict1[layer_name].get_weights()
	if w2_train is None: w2_train = [1.]*len(x2_train)
	if shuffle:
		m1_index = range(len(x1_train))
		for n in range(10): np.random.shuffle(m1_index)
		x1_train = x1_train[m1_index]
		y1_train = y1_train[m1_index]
		m2_index = range(len(x2_train))
		for n in range(10): np.random.shuffle(m2_index)
		x2_train = x2_train[m2_index]
		y2_train = y2_train[m2_index]
		w2_train = w2_train[m2_index]
	
	epoch = 1
	index = 0
	best_iter = 0
	best_auc = np.asarray([-np.inf, -np.inf])
	best_auc_m1 = np.asarray([-np.inf, -np.inf])
	best_auc_m2 = np.asarray([-np.inf, -np.inf])
	best_model = np.asarray([1,1])
	
	patience = 200
	patience_increase = 2
	improvement_threshold = 1.01
	validation_frequency = 1
	done_looping = False
	
	while (epoch < n_epoch) and (not done_looping):
		iter = index
		m1_idx = construct_balanced_minibatch(index, r10_index, r11_index, batch_size1)
		if (not fix_shared_layer) and enable_sharing:
			model1.train_on_batch(
				x1_train[m1_idx,:], 
				y1_train[m1_idx])
		if enable_sharing:
			for layer_name in shared_layer:
				layer_dict2[layer_name].set_weights(layer_dict1[layer_name].get_weights())
		m2_idx = construct_balanced_minibatch(index, r20_index, r21_index, batch_size2)
		model2.train_on_batch(
			x2_train[m2_idx,:], 
			y2_train[m2_idx],
			sample_weight=w2_train[m2_idx])
		if fix_shared_layer: 
			#layer_dict2['fx'].set_weights([fx_W, fx_b]) 
			for layer_name in shared_layer:
				layer_dict2[layer_name].set_weights(shared_layer[layer_name])
		if enable_sharing:
			for layer_name in shared_layer:
				layer_dict1[layer_name].set_weights(layer_dict2[layer_name].get_weights())
		
		if iter % 100 == 0:
				print('training @ iter = ' + str(iter))
		if not (iter+1) % validation_frequency: 
			this_y_val_hat11 = model1.predict(x1_val)
			this_y_val_hat12 = model1.predict(x2_val)
			this_val_eval11 = model1.evaluate(x1_val, y1_val)
			this_val_eval12 = model1.evaluate(x2_val, y2_val)
			this_y_val_auc11 = metrics.roc_auc_score((y1_val+1)/2., (this_y_val_hat11+1)/2.)
			this_y_val_auc12 = metrics.roc_auc_score((y2_val+1)/2., (this_y_val_hat12+1)/2.)
			auc1_list = np.asarray([this_y_val_auc11, this_y_val_auc12])
			
			this_y_val_hat21 = model2.predict(x1_val)
			this_y_val_hat22 = model2.predict(x2_val)
			this_val_eval21 = model2.evaluate(x1_val, y1_val)
			this_val_eval22 = model2.evaluate(x2_val, y2_val)
			this_y_val_auc21 = metrics.roc_auc_score((y1_val+1)/2., (this_y_val_hat21+1)/2.)
			this_y_val_auc22 = metrics.roc_auc_score((y2_val+1)/2., (this_y_val_hat22+1)/2.)
			auc2_list = np.asarray([this_y_val_auc21, this_y_val_auc22])
			
			this_auc = np.max([auc1_list, auc2_list], axis=0)
	
			print('epoch %i, minibatch %i/%i, validation auc [%.3f, %.3f]'%
				(epoch, (index + 1)%n_train_batches, n_train_batches, this_auc[0]*100, this_auc[1]*100)) 
			print(' auc1=[%.3f, %.3f], auc2=[%.3f, %.3f]'%
				(auc1_list[0]*100, auc1_list[1]*100, 
				auc2_list[0]*100, auc2_list[1]*100))
			
			if all(auc1_list>best_auc_m1):
				best_auc_m1 = np.copy(auc1_list)
				model1.save(save_model1_name)
			if all(auc2_list>best_auc_m2):
				best_auc_m2 = np.copy(auc2_list)
				model2.save(save_model2_name)
			
			if any(this_auc > best_auc):
				if any(this_auc > best_auc * improvement_threshold):
					patience = max(patience, iter*patience_increase)
				best_auc[this_auc>best_auc] = this_auc[this_auc>best_auc]
				best_iter = iter
				if any(best_auc == auc1_list):
					model1.save(save_model1_name)
					best_model[best_auc == auc1_list]=1
					best_auc_m1 = np.copy(auc1_list)
				if any(best_auc == auc2_list):
					model2.save(save_model2_name)
					best_model[best_auc == auc2_list]=2
					best_auc_m2 = np.copy(auc2_list)
				
				print "-"*50
				print epoch
				print ("     epoch %i, minibatch %i/%i") % (epoch, (index+1)%n_train_batches, n_train_batches)
				print "best auc = "+str(best_auc)
				print "-"*50
		if patience <= iter:
			done_looping=True
			break
			
		index += 1
		if index % n_train_batches ==0:
			epoch += 1
	
	model1.load_weights(save_model1_name)
	model2.load_weights(save_model2_name)
	
	print "Before ensemble: Best auc=%s @ of model=%s"%(str(best_auc*100), str(best_model))
	if ensemble:
		this_y_val_hat31 = (model1.predict(x1_val) + model2.predict(x1_val))/2
		this_y_val_hat32 = (model1.predict(x2_val) + model2.predict(x2_val))/2
		this_y_val_auc31 = metrics.roc_auc_score((y1_val+1)/2., (this_y_val_hat31+1)/2.)
		this_y_val_auc32 = metrics.roc_auc_score((y2_val+1)/2., (this_y_val_hat32+1)/2.)
		auc3_list = np.asarray([this_y_val_auc31, this_y_val_auc32])
		print auc3_list
		if any(auc3_list>best_auc): 
			best_model[auc3_list>best_auc]=3
			best_auc[auc3_list>best_auc]=auc3_list[auc3_list>best_auc]
			print "After ensemble: Best auc=%s @ of model=%s"%(str(best_auc*100), str(best_model))
	
	
	print "Optimization complete."
	print "Best individual model, m1=%s, m2=%s" % (str(best_auc_m1*100), str(best_auc_m2*100))
	print "Best auc=%s"%(str(best_auc*100))+ \
		' @ iter ' +str(best_iter) + \
		' of model='+str(best_model)
	
	return best_auc, best_model


def metrics_keras_model(model, x_val, y_val):
	""" given a keras model, return all possible metrics (ROC, PR, acc, )
	"""
	return


def train_keras_regression_model(model, x_train, y_train, x_val, y_val, sample_weight=None, 
		n_epoch=100, batch_size=64, save_model_name='my_model.h5', metric_to_use='loss'):	
	n_train_batches = len(x_train)//batch_size
	epoch = 1
	index = 0
	best_iter = 0
	best_auc = -np.inf if metric_to_use!='loss' else np.inf
	#metric_to_use = 'pr'
	
	patience = 1000
	patience_increase = 2
	improvement_threshold = 1.01 if metric_to_use!='loss' else 0.95
	validation_frequency = min(n_train_batches, patience //2)
	
	done_looping = False
	total_index = range(len(x_train))
	while (epoch < n_epoch) and (not done_looping):
		for index in range(n_train_batches):
			iter = (epoch - 1) * n_train_batches + index
			#sam_idx = total_index[(index*batch_size%len(x_train)):((index+1)*batch_size%len(x_train))]
			if (index*batch_size%len(x_train)) < ((index+1)*batch_size%len(x_train)):
				sam_idx = total_index[(index*batch_size%len(x_train)):((index+1)*batch_size%len(x_train))]
			else:
				sam_idx = total_index[(index*batch_size%len(x_train)):] + total_index[:((index+1)*batch_size%len(x_train))]
			if sample_weight is not None:
				model.train_on_batch(
					x_train[sam_idx,:], 
					y_train[sam_idx],
					sample_weight=sample_weight[sam_idx])
			else:
				model.train_on_batch(
					x_train[sam_idx,:], 
					y_train[sam_idx])
			if iter % 100 == 0:
				print('training @ iter = ' + str(iter))
			if (iter+1) % validation_frequency: 
				this_y_val_hat = model.predict(x_val)
				this_val_eval = model.evaluate(x_val, y_val)
				this_train_eval = model.evaluate(x_train, y_train)
				this_corr = scipy.stats.pearsonr(y_val, this_y_val_hat[:,0])[0]
				this_metrics = {'loss':this_val_eval,
					'corr':this_corr}
				this_auc = this_metrics[metric_to_use]
				print('epoch %i, minibatch %i/%i, train loss %.3f, validation loss %.3f, corr %.3f %%' %
					(epoch, index + 1, n_train_batches,
					this_train_eval,
					this_val_eval, this_corr))
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
		epoch += 1
	model.load_weights(save_model_name)
	print "Optimization complete. Best auc="+str(best_auc)+'@ iter ' +str(best_iter)
	return best_auc


def train_keras_weighted_model(model, x_train, y_train, x_val, y_val, 
		w_train, w_val, batch_size=64, n_epoch=100, 
		save_model_name='my_weighted_model.h5', metric_to_use='roc'):
	
	def eval_acc(y_true, y_pred, W2_val):
		corrected_predicted = 0
		for i in range(len(y_true)):
			if y_true[i]<=0.5 and y_pred[i]<=0:
				corrected_predicted+=W2_val[i]
			elif y_true[i]>0.5 and y_pred[i]>0:
				corrected_predicted+=W2_val[i]
		acc = float(corrected_predicted)/np.sum(W2_val)
		return acc
	
	r0_index = [i for i in range(len(y_train)) if y_train[i]==-1]
	r1_index = [i for i in range(len(y_train)) if y_train[i]==1]
	
	n_R1_batches = len(r1_index)//(batch_size/2)
	n_R0_batches = len(r0_index)//(batch_size/2)
	n_train_batches = max(n_R1_batches, n_R0_batches)
	epoch = 1
	index = 0
	best_iter = 0
	best_auc = -np.inf if metric_to_use!='loss' else np.inf
	#metric_to_use = 'pr'
	
	patience = 1000
	patience_increase = 2
	improvement_threshold = 1.01 if metric_to_use!='loss' else 0.95
	validation_frequency = min(n_train_batches, patience_increase //2)
	
	done_looping = False
	while (epoch < n_epoch) and (not done_looping):
		for index in range(n_train_batches):
			iter = (epoch - 1) * n_train_batches + index
			sam_idx = construct_balanced_minibatch(index, r0_index, r1_index, batch_size)
			
			model.train_on_batch(
				x_train[sam_idx,:], 
				y_train[sam_idx],
				sample_weight=w_train[sam_idx])
			
			if iter % 100 == 0:
				print('training @ iter = ' + str(iter))
			if not (iter+1) % validation_frequency: 
				this_y_val_hat = model.predict(x_val)
				this_val_eval = model.evaluate(x_val, y_val)
				this_train_eval = model.evaluate(x_train, y_train)
				this_roc = metrics.roc_auc_score((y_val+1)/2., (this_y_val_hat+1)/2., sample_weight=w_val)
				this_pr = metrics.average_precision_score((y_val+1)/2., (this_y_val_hat+1)/2., sample_weight=w_val)
				this_weighted_acc = eval_acc((y_val+1)/2., this_y_val_hat, w_val)
				this_metrics = {'acc':this_val_eval[1], 'loss':this_val_eval[0],
					'roc':this_roc, 'pr':this_pr, 'weighted_acc':this_weighted_acc}
				this_auc = this_metrics[metric_to_use]
				print('epoch %i, minibatch %i/%i, train loss %.3f, validation loss %.3f, corr %.3f %%' %
					(epoch, index + 1, n_train_batches,
					this_train_eval,
					this_val_eval, this_corr))
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
		epoch += 1
	model.load_weights(save_model_name)
	print "Optimization complete. Best auc="+str(best_auc)+'@ iter ' +str(best_iter)
	return best_auc


class boosting_ensemble_model(object):
	""" a class for adding new models and betas and do weighted sum prediction
	"""
	def __init__(self, model_list, beta_list, name_list=None):
		assert len(model_list)==len(beta_list)
		self.model_list = model_list
		self.beta_list = beta_list
		self.name_list = [None]*len(model_list) if name_list is None else name_list
	
	def add_model(self, model, beta, name=None):
		self.model_list.append(model)
		self.beta_list.append(beta)
		self.name_list.append(name)
	
	def delete_model(self, model_index):
		del self.model_list[model_index]
		del self.beta_list[model_index]
		del self.name_list[model_index]
	
	def ensemble_predict(self, X, n_out=1):
		pred = np.zeros((X.shape[0],n_out))
		for i in range(len(self.model_list)):
			pred += self.beta_list[i]*self.model_list[i].predict(X)
		return pred
	
	def print_formula(self):
		equ = []
		for i in range(len(self.model_list)):
			equ.append(str(self.beta_list[i]) + ' x ' + self.name_list[i])
		return ' + '.join(equ)