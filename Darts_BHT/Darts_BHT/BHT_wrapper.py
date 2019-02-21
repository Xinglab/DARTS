# -*- coding: UTF-8 -*-

import os
from .convert_rmats import read_rmats_counts, write_darts_counts_from_rmats
from . import pretty_writter
import logging

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
Darts_BHT_Rcpp = importr("Darts")

logger = logging.getLogger('Darts_BHT.bayes_infer')

def run_darts_BHT_norep(in_fn, out_fn, out_RData_fn, cutoff, rescale_meth, rho_fn, verbose, thread):
	frame = (in_fn,
		out_fn,
		out_RData_fn,
		rescale_meth,
		cutoff,
		rho_fn,
		verbose,
		thread)
	#print(frame)

	Darts_BHT_Rcpp.Darts(in_fn=in_fn,
		out_fn=out_fn,
		out_RData_fn=out_RData_fn,
		C=cutoff,
		rescale_meth=rescale_meth,
		rho_fn=rho_fn,
		verbose=verbose,
		thread=thread)
	return 0


def run_darts_BHT_rep(in_fn, out_fn, out_RData_fn, cutoff, rescale_meth, rho_fn, estim_gVar, is_paired, pooled, verbose, thread):
	frame = (in_fn,
		out_fn,
		out_RData_fn,
		rescale_meth,
		cutoff,
		rho_fn,
		estim_gVar,
		is_paired,
		pooled,
		verbose,
		thread)
	#print(frame)

	Darts_BHT_Rcpp.Darts_replicate(in_fn=in_fn,
		out_fn=out_fn,
		out_RData_fn=out_RData_fn,
		rescale_meth=rescale_meth,
		C=cutoff,
		rho_fn=rho_fn,
		estim_groupVar_prior=estim_gVar,
		is_paired=is_paired,
		pooling=pooled,
		verbose=verbose,
		thread=thread)
	return 0


def validate_count_file(fp, replicate_model):
	has_replicates = True
	with open(fp, 'r') as f:
		firstline = True
		for line in f:
			ele = line.strip().split()
			if firstline:
				header = {ele[i]:i for i in range(len(ele))}
				if 'IJC_SAMPLE_1' in ele and 'IJC_SAMPLE_2' in ele:
					count_names = ['IJC_SAMPLE_1', 'IJC_SAMPLE_2', 'SJC_SAMPLE_1', 'SJC_SAMPLE_2']
				else:
					count_names = ['I1', 'I2', 'S1', 'S2']
				firstline = False
				continue
			this_data = [ele[header[x]].split(',') for x in count_names]
			this_data_len = [x for x in map(len, this_data)]
			if not all([x>1 for x in this_data_len]):
				has_replicates = False
			if replicate_model=="paired" and len(set(this_data_len))!=1:
				logger.info('detected un-paired data at line %s'%line)
	return has_replicates


def parser(args):
	is_paired = args.replicate_model == 'paired'
	pooled = args.replicate_model == 'pooled'
	event_type = args.event_type
	rescale_meth_map = {'gaussian_mixture':1, 'bias_estimates':2}
	args.rescale_method = rescale_meth_map[args.rescale_method]
	if not os.path.isdir(args.outdir):
		os.makedirs(args.outdir)
	validated_count_fp = os.path.join(args.outdir, "{}.input.txt".format(args.event_type) )

	# parsing counts file
	if args.rmats_count_fp:
		if not os.path.isfile(args.rmats_count_fp):
			raise Exception("input count file '%s' is not found"%in_fn)

		logger.info('Coverting rMATS count to Darts format')
		exon_dict, _ = read_rmats_counts(count_fp=args.rmats_count_fp, annot_fp=args.annot, event_type=args.event_type)
		write_darts_counts_from_rmats(exon_dict, fn=validated_count_fp)
	elif args.darts_count_fp:
		validated_count_fp = args.darts_count_fp

	logger.info('input count={}'.format(validated_count_fp))
	logger.info('output dir={}'.format(args.outdir))

	# auto-detect if has replicates
	has_replicates = validate_count_file(validated_count_fp, args.replicate_model)
	
	if not has_replicates and (args.estim_gVar or args.replicate_model!='none'):
		logger.info('detected input file has no replicates; your "replicate_model" or "estim_gVar" options will be ignored')
	if has_replicates and args.replicate_model=="none":
		args.replicate_model='unpaired'

	# call Rcpp code
	prior_suffix = 'info' if args.prior else 'flat'
	out_fn = os.path.join(args.outdir, "{}.darts_bht.{}.txt".format(args.event_type, prior_suffix))
	out_RData_fn = os.path.join(args.outdir, "{}.darts_bht.{}.RData".format(args.event_type, prior_suffix))
	args.prior = ro.NA_Integer if not args.prior else args.prior
	if has_replicates:
		logger.info('using replicate model, mode "{}"'.format( args.replicate_model ) )
		run_darts_BHT_rep(
			in_fn=validated_count_fp,
			out_fn=out_fn,
			out_RData_fn=out_RData_fn,
			cutoff=args.cutoff,
			rescale_meth=args.rescale_method,
			rho_fn=args.prior,
			estim_gVar=args.estim_gVar,
			is_paired=is_paired,
			pooled=pooled,
			verbose=args.verbose,
			thread=args.nthread
			)
	else:
		logger.info('using no-replicate model')
		run_darts_BHT_norep(
			in_fn=validated_count_fp, 
			out_fn=out_fn,
			out_RData_fn=out_RData_fn,
			cutoff=args.cutoff, 
			rescale_meth=args.rescale_method, 
			rho_fn=args.prior, 
			verbose=args.verbose, 
			thread=args.nthread)

	## add module beatify results in XLSX format by integrating annotation file
	pretty_writter.write_xlsx(out_fn, args.annot, args.event_type)
	return