"""
# sklearn-style neural net in keras

# Zijun Zhang, 3/9/2017

# revised 6.6.2017: add batch normalization options
"""

from keras.utils import np_utils
from keras.models import Sequential, Model
from keras.layers import Dense, Dropout, Activation, Input
from keras.layers.normalization import BatchNormalization
from keras.optimizers import RMSprop, SGD
from keras.constraints import maxnorm
from keras.callbacks import ModelCheckpoint, EarlyStopping
from keras.callbacks import Callback

import numpy as np
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.utils.validation import check_X_y, check_array, check_is_fitted
#from sklearn.utils.multiclass import unique_labels
#from sklearn.metrics import euclidean_distances
from .utils import *

class KerasNeuralNetClassifier(BaseEstimator, ClassifierMixin):
	"""
	"""
	def __init__(self, n_in, n_out, batch_size=64, n_hidden=[16], init_func='glorot_normal',
	            dropout=[0, 0.5], activation='relu', optimizer='RMSprop', loss='categorical_crossentropy',
				n_epoches=100, fit_method='balanced', output_layer='softmax',
				batch_norm=True,
				pred_transform=None, label_transform=None):
		
		assert len(n_hidden)+1 == len(dropout)
		input = Input(shape=(n_in,), dtype='float32', name='input')
		firstlayer = True
		if dropout[0]>0:
			h = Dropout(dropout[0])(input)
			firstlayer = False
		for i in range(len(n_hidden)):
			if firstlayer:
				h = Dense(n_hidden[i], activation=None, init=init_func)(input)
				if batch_norm:
					h = BatchNormalization()(h)
				h = Activation(activation)(h)
				firstlayer=False
			else:
				h = Dense(n_hidden[i], activation=None, init=init_func)(h)
				if batch_norm:
					h = BatchNormalization()(h)
				h = Activation(activation)(h)
			if dropout[i+1]>0: 
				h = Dropout(dropout[i+1])(h)
		
		output = Dense(n_out, activation=output_layer, name='output')(h)
		model = Model(input=input, output=output)
		model.compile(optimizer=optimizer, loss=loss, metrics=['accuracy'])
		self.model = model
		self.batch_size = batch_size
		self.multi_output = True if n_out>1 else False
		self.n_in = n_in
		self.n_out = n_out
		self.n_hidden = n_hidden
		self.init_func = init_func
		self.dropout = dropout
		self.activation = activation
		self.optimizer = optimizer
		self.loss = loss
		self.n_epoches = n_epoches
		self.fit_method = fit_method
		self.output_layer = output_layer
		self.pred_transform = pred_transform
		self.label_transform = label_transform
	
	def __str__(self):
		self.model.summary()
		return ''
	
	def _fit_with_val(self, X, y):
		#X, y = check_X_y(X, y, self.multi_output)
		checkpointer = ModelCheckpoint(filepath= ".keras_model.h5", verbose=1, save_best_only=True, monitor='val_loss')
		earlystopper = EarlyStopping(monitor='val_loss', patience=100, verbose=1)
		history = self.model.fit(
			X, y, 
			batch_size=self.batch_size, 
			nb_epoch=self.n_epoches, 
			verbose=0,
			shuffle=True,
			callbacks=[checkpointer, earlystopper],
			validation_split=0.1
		)
		self.model.load_weights(".keras_model.h5")
		return self
	
	def _overfit(self, X, y):
		self.model.fit(
			X, y, 
			batch_size=self.batch_size, 
			nb_epoch=self.n_epoches, 
			verbose=0,
			shuffle=True
		)
		return self
		
	def _balanced_fit(self, X, y, nb_fold=5, metric_to_use='roc', minibatch_pos_prop=0.5, random_seed=123456):
		X_val, Y_val, X_train, Y_train = split_data_fold(X, y, nb_fold, 0, random_seed=random_seed)
		best_auc = train_keras_balanced_model(
			self.model, 
			X_train, Y_train, 
			X_val, Y_val, 
			metric_to_use=metric_to_use,
			batch_size=self.batch_size,
			n_epoch=self.n_epoches,
			minibatch_pos_prop=minibatch_pos_prop,
			save_model_name='.keras_model.h5'
		)
		return best_auc
	
	def fit(self, X, y, **kwargs):
		if self.label_transform is not None:
			y = self.label_transform(y)
		if self.fit_method=='overfit':
			res = self._overfit(X,y)
		elif self.fit_method=='keras_fit':
			res = self._fit_with_val(X,y)
		elif self.fit_method=='balanced':
			res = self._balanced_fit(X,y,**kwargs)
		return res
	
	def predict(self, X):
		X = check_array(X)
		if self.pred_transform is None:
			y_pred = self.model.predict(X)
		else:
			y_pred = self.pred_transform(self.model.predict(X))
		return y_pred
			