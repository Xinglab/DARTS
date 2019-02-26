#!/usr/bin/env python

from setuptools import setup


def main():
	setup(name='Darts_DNN',
		  
		  version='0.1.0',
		  
		  description='Deep-learning Augmented RNA-seq analyis of Transcript Splicing',
		  
		  author='Zijun Zhang',
		  
		  author_email='zj.z@ucla.edu',
		  
		  url='https://github.com/Xinglab/DARTS',
		  
		  packages=['Darts_DNN'],
		  
		  scripts=['bin/Darts_DNN']
		 )
	return

if __name__ == '__main__':
	main()