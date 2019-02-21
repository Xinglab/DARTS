'''use pandas to store Darts_BHT output
from text to xlsx
ZZJ
2.20.2019
'''

import os
import pandas as pd
from openpyxl import load_workbook
import logging

logger = logging.getLogger('Darts_BHT.bayes_infer')

rmats_type = {
	"SE": ["chr", "strand", "upstreamEE", "exonStart_0base", "exonEnd", "downstreamES"],
	"A3SS": ["chr", "strand", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"],
	"A5SS": ["chr", "strand", "longExonStart_0base", "longExonEnd",  "shortES", "shortEE", "flankingES", "flankingEE"],
	"RI": ["chr", "strand", "riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
	'MXE': ["chr", "strand", "1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"]
}


def write_xlsx(res_fp, annot_fp, event_type):
	xlsx_fp = os.path.join(os.path.dirname(res_fp), 'Darts_BHT.result.xlsx')
	## read in txt file
	res_df = pd.read_table(res_fp, index_col=0)
	annot_df = pd.read_table(annot_fp, index_col=0)
	newID = [':'.join([str(x) for x in annot_df.loc[i, rmats_type[event_type] ]]) for i in range(annot_df.shape[0]) ]
	annot_df.index = newID

	## process xlsx
	xlsx_df = pd.merge(res_df, annot_df, how='left', left_index=True, right_index=True)
	# DO SOME OTHER FANCY STUFF HERE
	# LIKE INSERT HYPERLINKS TO EXONS

	## write excel file
	if os.path.isfile(xlsx_fp):
		book = load_workbook(xlsx_fp)
		add_sheet = True
	else:
		add_sheet = False
	prior_type = 'flat' if 'flat' in os.path.basename(res_fp) else 'info'
	writer = pd.ExcelWriter(xlsx_fp, engine='openpyxl')
	if add_sheet:
		writer.book = book
	xlsx_df.to_excel(writer, sheet_name="{}-{}".format(event_type, prior_type))
	writer.save()
	writer.close()