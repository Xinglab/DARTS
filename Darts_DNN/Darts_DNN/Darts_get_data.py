# -*- coding: UTF-8 -*-

"""
Darts_DNN - get_data

Implements an internal downloading module for 
getting data from internet. Data includes training
data, cis feature files, etc.

This module depends on the url and md5sum stored in ``resources/download.yaml``


"""

import os
import logging

from . import download
from . import config

logger = logging.getLogger('Darts_DNN.get_data')


def parser(args):
	# parse options
	data = list(set(args.data))
	event_type = args.event_type
	outdir = args.out_dir
	data_config = config.DOWNLOAD_CONFIG[args.sp]

	# perform downloading
	for d in data:
		if not event_type in data_config[d]:
			logger.info("data '%s' is not available in event type '%s' for version '%s'"%(d, event_type, config.CURRENT_VERSION) )
			continue
		url_list = data_config[d][event_type]['url']
		md5_list = data_config[d][event_type]['md5sum']
		output_list = data_config[d][event_type]['output']
		url_list = url_list if isinstance(url_list, list) else [url_list]
		md5_list = md5_list if isinstance(md5_list, list) else [md5_list]
		output_list = output_list if isinstance(output_list, list) else [output_list]
		for url, md5, outfn in zip(url_list, md5_list, output_list):
			md5 = None if md5=='None' else md5
			if outdir or outfn == 'None':
				basename = os.path.basename(url)
				outfn = os.path.join(outdir, basename) if outdir else os.path.join('.', basename)
			if outfn.startswith('~/'):
				outfn = os.path.join(os.path.expanduser('~'), outfn.lstrip('~/'))
			if not os.path.isdir(os.path.dirname(outfn)):
				os.makedirs(os.path.dirname(outfn))
			logger.info("dowloading '%s' to '%s'" % (url, outfn) )
			download.download_with_resume(url, outfn, md5)
	return
