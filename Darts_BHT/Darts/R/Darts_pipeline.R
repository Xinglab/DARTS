#!/usr/bin/env Rscript

# pipeline utilities for Darts_BHT
# Zijun Zhang
# 5.14.2017
# revised 7.1.2017: modified `Spider` to accept informative prior
# revised 7.2.2017: added options and a function to rescale informative prior
# revised 7.14.2017: changed name to `Darts` and built into a package
# revised 7.17.2017: put this wrapper within package
# revised 8.20.2017: added random_state for reproducibility
# revised 3.19.2018: moved rescaling-related functions to `rescale_by_bias.R`
# revised 3.19.2018: use new bias-estimates to rescale


match_input_column_names = function(df)
{
	rmats_names = c('ID', 'IJC_SAMPLE_1', 'SJC_SAMPLE_1', 'IJC_SAMPLE_2', 'SJC_SAMPLE_2', 'IncFormLen', 'SkipFormLen')
	darts_names = c('ID', 'I1', 'S1', 'I2', 'S2', 'inc_len', 'skp_len')
	name_matching = data.frame(rmats=rmats_names, darts=darts_names)
	idx = match(name_matching$rmats, colnames(df), nomatch=0)
	df = df[,idx]
	idx2 = match(colnames(df), name_matching$rmats)
	colnames(df) = name_matching$darts[idx2]
	return(df)
}


check_data_sanity = function(data)
{
	is_healthy=1
	if(data$I1==0 && data$S1==0) is_healthy=0
	if(data$I2==0 && data$S2==0) is_healthy=0
	#if( min(data$I1, data$S1)<=2 ) is_healthy=0
	#if( min(data$I2, data$S2)<=2 ) is_healthy=0
	#if( sum(data$I1, data$S1)<=20 ) is_healthy=0
	#if( sum(data$I2, data$S2)<=20 ) is_healthy=0
	#if(data$I1==0 && data$I2==0) is_healthy=0
	#if(data$S1==0 && data$S2==0) is_healthy=0
	return(is_healthy)
}

check_data_sanity_replicate = function(data)
{
	is_healthy=1
	n_rep = nrow(data)/2
	if(sum(data$Inc==0 & data$Skp==0)>0) is_healthy=0
	if(all(data$Skp==0)) is_healthy=0
	if(all(data$Inc==0)) is_healthy=0
	#if(all(data$Skp[data$j==1]==0)) is_healthy=0
	#if(all(data$Skp[data$j==2]==0)) is_healthy=0
	#if(all(data$Inc[data$j==1]==0)) is_healthy=0
	#if(all(data$Inc[data$j==2]==0)) is_healthy=0
	return(is_healthy)
}

change_cutoff = function(in_fn, rdata_fn, C)
{
	data=read.table(in_fn, sep='\t', header=T)
	load(rdata_fn)
	delta_quantiles = seq(-0.99, 0.99, 0.01)
	right = which(abs(delta_quantiles-abs(C))<0.0001)
	left = which(abs(delta_quantiles+abs(C))<0.0001)
	post_pr = 1-res_dict[,right] + res_dict[,left]
	data$post_pr = post_pr
	return(data)
}


