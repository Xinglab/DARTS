#!/usr/bin/env python

from setuptools import setup


def main():
	setup(
		name='Darts_DNN',
		  
		version='0.1.0',
		  
		description='Deep-learning Augmented RNA-seq analysis of Transcript Splicing',

		author='Zijun Zhang',

		author_email='zj.z@ucla.edu',

		url='https://github.com/Xinglab/DARTS',

		packages=['Darts_DNN', 'Darts_DNN.resources'],

		scripts=['bin/Darts_DNN'],

		include_package_data=True,

		package_data={'Darts_DNN.resources':[
				'architecture.yaml', 
		  		'cis_features.yaml', 
				'download.yaml',
				'human_t2g.txt',
				'rbp_gene_list.txt',
				'trained_model_param.yaml',
				'trans_features.yaml']},
		install_requires=[]#'keras', 
				#'numpy',
				#'pyyaml',
				#'h5py',
				#'scikit-learn',
				#'scipy',
				#'tqdm>=4.14',
				#'pandas>=0.21.0',
				#'theano']
		 )
	return

if __name__ == '__main__':
	main()
