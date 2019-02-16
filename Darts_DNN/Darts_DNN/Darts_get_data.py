# -*- coding: UTF-8 -*-

"""
Darts_DNN - get_data

date: Feb. 15, 2019
author: ZZJ

Implements an internal downloading module for 
getting data from internet. Data includes training
data, cis feature files, etc.
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

	# perform downloading
	for d in data:
		if not event_type in config.CURRENT_AVAILABLE_DATA[d]:
			logger.info("data '%s' is not available in event type '%s' for version '%s'"%(d, event_type, config.CURRENT_VERSION) )
			continue
		url = config.CURRENT_AVAILABLE_DATA[d][event_type]['url']
		md5 = config.CURRENT_AVAILABLE_DATA[d][event_type]['md5sum']
		md5 = None if md5=='None' else md5
		basename = os.path.basename(url)
		outfn = os.path.join(outdir, basename)
		logger.info("dowloading '%s' to '%s'" % (url, outfn) )
		download.download_with_resume(url, outfn, md5)
	return
