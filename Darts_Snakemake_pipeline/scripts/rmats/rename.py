# make the count table for all Control samples
# between var1 and var2
# Zijun Zhang
# 6.15.2017
# revised 8.12.2017: changed to Kallisto results

import sys
import os
from collections import defaultdict
from datetime import datetime
import random

def read_rmats_counts(dir):
	exon_dict = defaultdict(dict)
	exon_id_to_eid = {}
	eid_list = []
	with open(os.path.join(dir, 'fromGTF.SE.txt'), 'r') as f:
		firstline = True
		for line in f:
			ele = line.rstrip().split()
			if firstline:
				header = {ele[x]:x for x in range(len(ele))}
				firstline=False
				continue
			chr = ele[header['chr']]
			strand = ele[header['strand']]
			start = ele[header['exonStart_0base']]
			end = ele[header['exonEnd']]
			up_end = ele[header['upstreamEE']]
			dw_start = ele[header['downstreamES']]
			eid = ':'.join([chr, strand, up_end, start, end, dw_start])
			id = ele[header['ID']]
			exon_id_to_eid[id] = eid
			eid_list.append(eid)
	with open(os.path.join(dir, 'JC.raw.input.SE.txt'), 'r') as f:
		firstline=True
		for line in f:
			ele=line.rstrip().split('\t')
			if firstline:
				header = {ele[x]:x for x in range(len(ele))}
				firstline=False
				continue
			id = ele[header['ID']]
			I1 = ele[header['IJC_SAMPLE_1']]
			I2 = ele[header['IJC_SAMPLE_2']]
			S1 = ele[header['SJC_SAMPLE_1']]
			S2 = ele[header['SJC_SAMPLE_2']]
			inc_len = int(ele[header['IncFormLen']])
			skp_len = int(ele[header['SkipFormLen']])
			exon_dict[exon_id_to_eid[id]]['I1'] = I1
			exon_dict[exon_id_to_eid[id]]['I2'] = I2
			exon_dict[exon_id_to_eid[id]]['S1'] = S1
			exon_dict[exon_id_to_eid[id]]['S2'] = S2
			exon_dict[exon_id_to_eid[id]]['inc_len'] = inc_len
			exon_dict[exon_id_to_eid[id]]['skp_len'] = skp_len
	return exon_dict, eid_list


def write_count_file(exon_dict, eid_list):
	#fout.write('\t'.join(['ID', 'I1', 'S1', 'I2', 'S2', 'inc_len', 'skp_len'])+'\n')
	print('\t'.join(['ID', 'I1', 'S1', 'I2', 'S2', 'inc_len', 'skp_len']))
	#for eid in exon_dict:
	for eid in eid_list:
		if not eid in exon_dict:
			continue
		I1 = exon_dict[eid]['I1']
		I2 = exon_dict[eid]['I2']
		S1 = exon_dict[eid]['S1']
		S2 = exon_dict[eid]['S2']
		inc_len = exon_dict[eid]['inc_len']
		skp_len = exon_dict[eid]['skp_len']
		#fout.write('\t'.join([eid, str(I1), str(S1), str(I2), str(S2), str(inc_len), str(skp_len)])+'\n')
		print('\t'.join([eid, str(I1), str(S1), str(I2), str(S2), str(inc_len), str(skp_len)]))
	return


if __name__ == '__main__':
	fp = sys.argv[1]
	dir = os.path.dirname(fp)
	exon_dict, eid_list = read_rmats_counts(dir)
	write_count_file(exon_dict, eid_list)