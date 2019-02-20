# -*- coding: UTF-8 -*-

'''Utils for connecting rMATS to Darts
'''

from collections import defaultdict

rmats_type = {
	"SE": ["chr", "strand", "upstreamEE", "exonStart_0base", "exonEnd", "downstreamES"],
	"A3SS": ["chr", "strand", "longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"],
	"A5SS": ["chr", "strand", "longExonStart_0base", "longExonEnd",  "shortES", "shortEE", "flankingES", "flankingEE"],
	"RI": ["chr", "strand", "riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
	'MXE': ["chr", "strand", "1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"]
}


def read_rmats_counts(count_fp, annot_fp, event_type='SE'):
	exon_dict = defaultdict(lambda: defaultdict(int))
	exon_id_to_eid = {}
	with open(os.path.join(dir, 'fromGTF.%s.txt'%event_type), 'r') as f:
		firstline = True
		for line in f:
			ele = line.rstrip().split()
			if firstline:
				header = {ele[x]:x for x in range(len(ele))}
				firstline=False
				continue
			eid = ':'.join([ ele[header[x]] for x in rmats_type[event_type] ])
			id = ele[header['ID']]
			exon_id_to_eid[id] = eid
	has_replicates = False
	with open(os.path.join(dir, 'JC.raw.input.%s.txt'%event_type), 'r') as f:
		firstline=True
		for line in f:
			ele=line.rstrip().split('\t')
			if firstline:
				header = {ele[x]:x for x in range(len(ele))}
				firstline=False
				continue
			id = ele[header['ID']]
			Inc1 = [int(x) for x in ele[header['IJC_SAMPLE_1']].split(',')]
			Inc2 = [int(x) for x in ele[header['IJC_SAMPLE_2']].split(',')]
			Skp1 = [int(x) for x in ele[header['SJC_SAMPLE_1']].split(',')]
			Skp2 = [int(x) for x in ele[header['SJC_SAMPLE_2']].split(',')]
			assert len(Inc1)==len(Skp1)
			assert len(Inc2)==len(Skp2)
			if len(Inc1)>1 and len(Inc2)>1:
				has_replicates = True
			inc_len = int(ele[header['IncFormLen']])
			skp_len = int(ele[header['SkipFormLen']])
			exon_dict[exon_id_to_eid[id]]['Inc1'] = Inc1
			exon_dict[exon_id_to_eid[id]]['Inc2'] = Inc2
			exon_dict[exon_id_to_eid[id]]['Skp1'] = Skp1
			exon_dict[exon_id_to_eid[id]]['Skp2'] = Skp2
			exon_dict[exon_id_to_eid[id]]['inc_len'] = inc_len
			exon_dict[exon_id_to_eid[id]]['skp_len'] = skp_len
	return exon_dict, has_replicates


def write_darts_counts_from_rmats(exon_dict, fn):
	with open(fn, 'w') as fout:
			fout.write('\t'.join(['ID', 'I1', 'S1', 'I2', 'S2', 'inc_len', 'skp_len'])+'\n')
			for eid in exon_dict:
					I1 = exon_dict[eid]['Inc1']
					S1 = exon_dict[eid]['Skp1']
					I2 = exon_dict[eid]['Inc2']
					S2 = exon_dict[eid]['Skp2']
					inc_len = target_exon_dict[eid]['inc_len']
					skp_len = target_exon_dict[eid]['skp_len']
					fout.write('\t'.join([eid, str(I1), str(S1), str(I2), str(S2), str(inc_len), str(skp_len)])+'\n')
	return