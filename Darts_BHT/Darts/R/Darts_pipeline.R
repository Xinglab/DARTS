#!/usr/bin/env Rscript

# the equivalent of SP_pipeline.py in R
# because Rcpp seems to run faster in R compared to rpy2
# Zijun Zhang
# 5.14.2017
# revised 7.1.2017: modified `Spider` to accept informative prior
# revised 7.2.2017: added options and a function to rescale informative prior
# revised 7.14.2017: changed name to `Darts` and built into a package
# revised 7.17.2017: put this wrapper within package
# revised 8.20.2017: added random_state for reproducibility

#library('Darts')

Darts = function(in_fn, out_dir, C=0.05, rho_fn=NA, verbose=1, random_state=777)
{
	out_fn = file.path(out_dir, ifelse(is.na(rho_fn), 'Sp_out.RData', 'Sp_out.prior.RData'))
	if(class(in_fn)=="character") {
		data=read.table(in_fn, sep='\t', header=T)
	} else if(class(in_fn)=="data.frame")
	{
		data = in_fn
	}
	if("IJC_SAMPLE_1" %in% colnames(data)) {
		data = match_input_column_names(data)
	}
	if(class(rho_fn)=="character")
	{
		rho.df = read.table(rho_fn, header=T)
		rho.df$rho = rescale_binary_prediction(rho.df$Y_pred, mu0.new=0.05, mu1.new=0.95, sd0.new=0.1, sd1.new=0.1, boundary=0.05)
		idx = match(data$ID, rho.df$ID)
		data$rho = rho.df$rho[idx]
	} else {
		if(! "rho" %in% colnames(data))
			data$rho = rep(0.5, nrow(data))
	}
	res = data
	res$mu.mle = NA
	res$delta.mle = NA
	res$post_pr = NA
	delta_quantiles = seq(-0.99, 0.99, 0.01)
	delta_quantiles = as.vector(delta_quantiles)
	res_dict = matrix(NA, nrow=nrow(data), ncol=length(delta_quantiles))
	rownames(res_dict) = data$exon_id
	
	pb = utils::txtProgressBar(min = 0, max = nrow(data), initial = 0) 
	for(i in 1:nrow(data))
	{
		if(verbose==2 && ! i%%200) cat(format(Sys.time(), "%c"), i,'/',nrow(data), '\n')
		if(verbose==1) utils::setTxtProgressBar(pb,i)
		id=data$ID[i]
		I1=data$I1[i]
		S1=data$S1[i]
		I2=data$I2[i]
		S2=data$S2[i]
		this_rho = data$rho[i]
		this_rho = ifelse(is.na(this_rho), 0.5, this_rho)
		inc_len = data$inc_len[i]
		skp_len = data$skp_len[i]
		this_data = list(I1=I1, S1=S1, I2=I2, S2=S2)
		is_healthy = check_data_sanity(this_data)
		if(! is_healthy) {
			next
		}
		set.seed(random_state)
		post = darts_bayesian_gibbs(rho=this_rho, data=this_data, 
			tau1=0.3, tau2=0.03,
			inc_eff_len=inc_len, skp_eff_len=skp_len)
		set.seed(NULL)
		post_cdf = sapply(delta_quantiles, function(x) mean(post[,2]<x))
		res_dict[i,] = post_cdf
		right = which(abs(delta_quantiles-abs(C))<0.0001)
		left = which(abs(delta_quantiles+abs(C))<0.0001)
		res$post_pr[i] = 1-post_cdf[right] + post_cdf[left]
		res$mu.mle[i] = I1/inc_len / (I1/inc_len + S1/skp_len)
		res$delta.mle[i] = I2/inc_len / (I2/inc_len + S2/skp_len) - I1/inc_len / (I1/inc_len + S1/skp_len)
	}
	close(pb)
	save(res_dict, file=out_fn)
	return(res)
}


Darts_replicate = function(in_fn, out_dir, C=0.05, rho_fn=NA, 
	estim_groupVar_prior=TRUE, is_paired=FALSE, pooling=FALSE, verbose=1, random_state=777)
{
	out_fn = file.path(out_dir, ifelse(is.na(rho_fn), 'Sp_out.RData', 'Sp_out.prior.RData'))
	if(class(in_fn)=="character") {
		data=read.table(in_fn, sep='\t', header=T, stringsAsFactors=F)
	} else if(class(in_fn)=="data.frame")
	{
		data = in_fn
	}
	if("IJC_SAMPLE_1" %in% colnames(data)) {
		data = match_input_column_names(data)
	}
	if(class(rho_fn)=="character")
	{
		rho.df = read.table(rho_fn, header=T)
		rho.df$rho = rescale_binary_prediction(rho.df$Y_pred, mu0.new=0.05, mu1.new=0.95, sd0.new=0.1, sd1.new=0.1, boundary=0.05)
		idx = match(data$ID, rho.df$ID)
		data$rho = rho.df$rho[idx]
	} else {
		if(! "rho" %in% colnames(data))
			data$rho = rep(0.5, nrow(data))
	}
	res = data
	res$mu.mle = NA
	res$delta.mle = NA
	res$post_pr = NA
	delta_quantiles = seq(-0.99, 0.99, 0.01)
	delta_quantiles = as.vector(delta_quantiles)
	res_dict = matrix(NA, nrow=nrow(data), ncol=length(delta_quantiles))
	rownames(res_dict) = data$exon_id
	
	if(estim_groupVar_prior) {
		prior_fit = estim_group_var.prior_fit(data)
		sgm_prior_fit = prior_fit$fit
		ggplot2::ggsave(file.path(out_dir,'groupVar.pdf'), plot=prior_fit$p, width=7, height=7)
	} else {
			sgm_prior_fit = NA
	}
	pb = utils::txtProgressBar(min = 0, max = nrow(data), initial = 0) 
	for(i in 1:nrow(data))
	{
		if(verbose==2 && ! i%%200) cat(format(Sys.time(), "%c"), i,'/',nrow(data), '\n')
		if(verbose==1) utils::setTxtProgressBar(pb,i)
		id=data$ID[i]
		I1=data$I1[i]; I1 = as.numeric(strsplit(as.character(I1),',')[[1]])
		S1=data$S1[i]; S1 = as.numeric(strsplit(as.character(S1),',')[[1]])
		I2=data$I2[i]; I2 = as.numeric(strsplit(as.character(I2),',')[[1]])
		S2=data$S2[i]; S2 = as.numeric(strsplit(as.character(S2),',')[[1]])
		this_rho = data$rho[i]
		this_rho = ifelse(is.na(this_rho), 0.5, this_rho)
		inc_len = data$inc_len[i]
		skp_len = data$skp_len[i]
		
		inc = c(I1, I2)
		skp = c(S1, S2)
		if(is_paired) {
			dat.k=rep(1:length(I1), 2)
		} else {
			dat.k=seq(1,length(inc))
		}
		dat.j = c(rep(1, length(I1)), rep(2, length(I2)))
		this_data = data.frame(k=dat.k, j=dat.j, Inc=inc, Skp=skp)
		set.seed(random_state)
		if(pooling) { 
			this_data = collapse_replicates(this_data)
			is_healthy = check_data_sanity(this_data)
			if(! is_healthy) next
			post = darts_bayesian_gibbs(rho=this_rho, data=this_data, 
				tau1=0.3, tau2=0.03,
				inc_eff_len=inc_len, skp_eff_len=skp_len)
			res$I1[i] = this_data$I1
			res$I2[i] = this_data$I2
			res$S1[i] = this_data$S1
			res$S2[i] = this_data$S2
		
		} else {
			is_healthy = check_data_sanity_replicate(this_data)
			if(! is_healthy) {
				next
			}
			post = rdarts_bayesian_gibbs(rho=this_rho, data=this_data, 
				tau1=0.3, tau2=0.03, sgm_prior_fit=sgm_prior_fit,
				inc_eff_len=inc_len, skp_eff_len=skp_len)
		}
		set.seed(NULL)
		post_cdf = sapply(delta_quantiles, function(x) mean(post[,2]<x))
		res_dict[i,] = post_cdf
		right = which(abs(delta_quantiles-abs(C))<0.0001)
		left = which(abs(delta_quantiles+abs(C))<0.0001)
		res$post_pr[i] = 1-post_cdf[right] + post_cdf[left]
		res$mu.mle[i] = sum(I1)/inc_len / (sum(I1)/inc_len + sum(S1)/skp_len)
		res$delta.mle[i] = sum(I2)/inc_len / (sum(I2)/inc_len + sum(S2)/skp_len) - sum(I1)/inc_len / (sum(I1)/inc_len + sum(S1)/skp_len)
	}
	close(pb)
	save(res_dict, file=out_fn)
	return(res)
}


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


rescale_binary_prediction = function(y_pred, 
	mu0.new=0.1, mu1.new=0.9, 
	sd0.new=NA, sd1.new=NA,
	boundary=0.1)
{
	#library("mixtools")
	rescale = function(OldValue, NewMax, NewMin, OldMax=1, OldMin=0)
	{
		OldRange = (OldMax - OldMin)  
		NewRange = (NewMax - NewMin)  
		NewValue = (((OldValue - OldMin) * NewRange) / OldRange) + NewMin
		return(NewValue)
	}
	mixmdl = mixtools::normalmixEM(y_pred, k=2)
	post = mixmdl$posterior
	if(mixmdl$mu[1]<mixmdl$mu[2]) {
		mu0 = mixmdl$mu[1]; sd0 = mixmdl$sigma[1]
		mu1 = mixmdl$mu[2]; sd1 = mixmdl$sigma[2]
	} else {
		mu0 = mixmdl$mu[2]; sd0 = mixmdl$sigma[2]
		mu1 = mixmdl$mu[1]; sd1 = mixmdl$sigma[1]
	}
	if(is.na(sd0.new)) sd0.new = sd0
	if(is.na(sd1.new)) sd1.new = sd1
	new_y_pred = c()
	for(i in 1:length(y_pred))
	{
		y = y_pred[i]
		q.0 = pnorm(y, mu0, sd0)
		q.1 = pnorm(y, mu1, sd1)
		q.0 = ifelse(q.0>0.99, 0.99, q.0)
		q.0 = ifelse(q.0<0.01, 0.01, q.0)
		q.1 = ifelse(q.1>0.99, 0.99, q.1)
		q.1 = ifelse(q.1<0.01, 0.01, q.1)
		score.0 = qnorm(q.0, mu0.new, sd0.new)
		score.1 = qnorm(q.1, mu1.new, sd1.new)
		if(mixmdl$mu[1]<mixmdl$mu[2])
		{
			new_pred = post[i,1]*score.0+post[i,2]*score.1
		} else {
			new_pred = post[i,2]*score.0+post[i,1]*score.1
		}
		#new_pred = ifelse(new_pred>0.99, 0.99, new_pred)
		#new_pred = ifelse(new_pred<0.01, 0.01, new_pred)
		new_y_pred = c(new_y_pred, new_pred)
	}
	new_y_pred = rescale(new_y_pred, 1-boundary, boundary, max(new_y_pred), min(new_y_pred))
	return(new_y_pred)
}

## if calling from Rscript with arguments, run command-line usage
# argv = commandArgs(trailingOnly=T)
# if(length(argv)>0)
# {
	# library("getopt")
	# spec <- matrix(c(
		# 'input'  , 'i', 1, "character", "input file in Spider format (required)",
		# 'rho'    , 'r', 2, "character", "prior file in Spider format (optional, if none then flat prior i.e. rho=0.5)",
		# 'out'    , 'o', 1, "character", "output directory (required)",
		# 'cutoff' , 'c', 2, "double",    "Cutoff of Posterior probability, must be in (0,1) (optional, default=0.05)",
		# 'rep'    , 'k', 2, "integer",   "run replicate [pair/unpair] model (optional, 0:no replicate, 1:unpaired, 2:paired, default=0)",
		# 'verbose', 'v', 2, "integer",   "verbose mode (optional, default=0)",
		# 'help'   , 'h', 0, "logical",   "print this help message and exit"
	# ),ncol=5,byrow=T)

	# opt = getopt(spec)
	# quit_session = function() { cat(paste(getopt(spec, usage=T),"\n"));q(); }
	# # assign default values
	# if ( is.null(opt$verbose ) ) { opt$verbose = FALSE }
	# if ( is.null(opt$cutoff ) ) { opt$cutoff = 0.05 }
	# if ( is.null(opt$rho ) ) { opt$rho = NA }
	# if ( is.null(opt$rep ) ) {opt$rep = 0}
	# ## check sanity of input arguments
	# if ( !is.null(opt$help) ) { quit_session() }
	# if( is.null(opt$input) || is.null(opt$out) || !file.exists(opt$input) ) {
		# cat("Either you didn't specify input/output filename or the input file doesn't exist.\n")
		# quit_session()
	# }
	# if( !is.na(opt$rho) && !file.exists(opt$rho)) {
		# cat("You specified a prior file:\n",opt$rho,",\n that does not exist.\n")
		# quit_session()
	# }
	
	# if( opt$cutoff<0 || opt$cutoff>1 ) {
		# cat("Cutoff value out of bounds.\n")
		# quit_session()
	# }
	
	# if( opt$rep<0 || opt$rep>2 ) {
		# cat("Invalid replicate model option.\n")
		# quit_session()
	# }
	
	# ## call stats model
	# if(is.na(opt$rho)) {
		# out_rdata_fn = file.path(opt$out, "Sp_out.RData")
		# out_table_fn = file.path(opt$out, "Sp_out.txt")
	# } else {
		# out_rdata_fn = file.path(opt$out, "Sp_out.prior.RData")
		# out_table_fn = file.path(opt$out, "Sp_out.prior.txt")
	# }
	# if( opt$k==0 ) {
		# res = Darts(opt$input, out_rdata_fn, C=opt$cutoff, rho_fn=opt$rho, verbose=opt$verbose)
	# } else
	# {
		# is_paired = ifelse(opt$k==2, TRUE, FALSE)
		# res = Darts_replicate(opt$input, out_rdata_fn, C=opt$cutoff, rho_fn=opt$rho, is_paired=is_paired, verbose=opt$verbose)
	# }
	# write.table(res, file=out_table_fn, quote=F, row.names=F, sep='\t')
# }