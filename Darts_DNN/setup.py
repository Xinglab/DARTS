#!/usr/bin/env python

from setuptools import setup


def main():
	setup(name='Darts_DNN',
		  
		  version='0.1.0',
		  
		  description='Deep-learning Augmented RNA-seq analyis of Transcript Splicing',
		  
		  author='Zijun Zhang',
		  
		  author_email='zj.z@ucla.edu',
		  
		  url='https://github.com/Xinglab/DARTS',
		  
		  packages=['Darts_DNN', 'Darts_DNN.resources'],
		  
		  scripts=['bin/Darts_DNN'],

		  include_package_data=True,

		  package_data={'Darts_DNN.resources':['rbp_max.tsv']}
		 )
	return

if __name__ == '__main__':
	main()