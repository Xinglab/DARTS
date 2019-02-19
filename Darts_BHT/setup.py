#!/usr/bin/env python

from setuptools import setup


def main():
	setup(
		name='Darts_BHT',
		  
		version='0.1.0',
		  
		description='Deep-learning Augmented RNA-seq analysis of Transcript Splicing',

		author='Zijun Zhang',

		author_email='zj.z@ucla.edu',

		url='https://github.com/Xinglab/DARTS',

		packages=['Darts_BHT'],

		scripts=['bin/Darts_BHT'],

		install_requires=['rpy2',
			'cython==0.27.0']
		 )

	from rpy2.robjects.packages import importr
	utils = importr('utils')
	utils.chooseCRANmirror(ind=1) 
	dependency_list = ['Rcpp', 'doSNOW', 'getopt', 'mixtools']
	#for package_name in dependency_list:
	#	utils.install_packages(package_name)
	utils.install_packages('Darts_BHT/Darts_0.1.0.tar.gz')
	return

if __name__ == '__main__':
	main()
