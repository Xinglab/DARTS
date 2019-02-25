# -*- coding: utf-8 -*-

'''Use pandas to store Darts_BHT output
from text to xlsx. If `openpyxl` is available, 
then append new sheets to existing xlsx file, making one "Darts_BHT.results.xlsx"
with multiple sheets, each sheet for a splicing event type;
otherwise, the existing xlsx file will be overwrite.
'''

import os
import pandas as pd
try:
	from openpyxl import load_workbook
	has_openpyxl = True
except:
	has_openpyxl = False

import logging

logger = logging.getLogger('Darts_BHT.result_writter')

rmats_type = {
	"SE": ["chr", "strand", "upstreamEE", "exonStart_0base", "exonEnd", "downstreamES"],
	"A3SS": ["chr", "strand", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"],
	"A5SS": ["chr", "strand", "longExonStart_0base", "longExonEnd",  "shortES", "shortEE", "flankingES", "flankingEE"],
	"RI": ["chr", "strand", "riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
	'MXE': ["chr", "strand", "1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"]
}


def write_xlsx(res_fp, annot_fp, event_type):
	xlsx_fp = os.path.join(os.path.dirname(res_fp), 'Darts_BHT.results.xlsx')
	## read in txt file
	logger.info('reading txt input')
	res_df = pd.read_csv(res_fp, index_col=0, sep='\t')
	res_df.rename({
		'inc_len':'IncFormLen', 
		'skp_len':'SkpFormLen',
		'rho': 'Prior',
		'psi1': 'PSI1',
		'psi2': 'PSI2',
		'mu.mle': 'AvgBasePSI',
		'delta.mle': 'IncLevelDiff',
		'post_pr': 'Posterior'}, axis='columns', inplace=True)
	annot_df = pd.read_csv(annot_fp, index_col=0, sep='\t')
	newID = [':'.join([str(x) for x in annot_df.loc[i, rmats_type[event_type] ]]) for i in range(annot_df.shape[0]) ]
	annot_df.index = newID

	## process xlsx
	xlsx_df = pd.merge(res_df, annot_df, how='left', left_index=True, right_index=True)
	# DO SOME OTHER FANCY STUFF HERE
	# LIKE INSERT HYPERLINKS TO EXONS

	xlsx_df = xlsx_df.sort_values(by='Posterior', ascending=False)

	## write excel file
	if os.path.isfile(xlsx_fp) and has_openpyxl:
		book = load_workbook(xlsx_fp)
		add_sheet = True
		logger.info('found previous results, appending a new sheet')
	else:
		add_sheet = False
		if os.path.isfile(xlsx_fp):
			logger.info('Overwrite Warning: found previous results, but openpyxl not available, will overwrite the file. Use "pip install openpyxl" to avoid this warning.')
	logger.info('writting xlsx')
	prior_type = 'flat' if 'flat' in os.path.basename(res_fp) else 'info'
	writer = pd.ExcelWriter(xlsx_fp, engine='openpyxl')
	if add_sheet:
		writer.book = book
	xlsx_df.to_excel(writer, sheet_name="{}-{}".format(event_type, prior_type))
	writer.save()
	writer.close()