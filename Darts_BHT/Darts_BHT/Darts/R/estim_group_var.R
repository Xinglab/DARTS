## estimate the within group variance by 
## fitting a loess regression model on read counts
## and imposing an inverse-gamma prior with pseudo-
## observation=1 (default).
# Zijun Zhang
# 7.14.2017

estim_group_var.prior_fit = function(data) 
{
	obs = matrix(NA, nrow=nrow(data), ncol=2)
	for(i in 1:nrow(data))
	{
		id=data$ID[i]
		I1=data$I1[i]; I1 = as.numeric(strsplit(as.character(I1),',')[[1]])
		S1=data$S1[i]; S1 = as.numeric(strsplit(as.character(S1),',')[[1]])
		I2=data$I2[i]; I2 = as.numeric(strsplit(as.character(I2),',')[[1]])
		S2=data$S2[i]; S2 = as.numeric(strsplit(as.character(S2),',')[[1]])
		inc_eff_len = data$inc_len[i]
		skp_eff_len = data$skp_len[i]
		
		inc = c(I1, I2)
		skp = c(S1, S2)
		dat.k=seq(1,length(inc))
	
		dat.j = c(rep(1, length(I1)), rep(2, length(I2)))
		this_data = data.frame(k=dat.k, j=dat.j, Inc=inc, Skp=skp)
		cond_idx = which(this_data$j==1)
		psi1 = this_data$Inc[cond_idx]/inc_eff_len / (this_data$Inc[cond_idx]/inc_eff_len + this_data$Skp[cond_idx]/skp_eff_len)
		psi2 = this_data$Inc[-cond_idx]/inc_eff_len / (this_data$Inc[-cond_idx]/inc_eff_len + this_data$Skp[-cond_idx]/skp_eff_len)
		
		sigma2 = (var(psi1) + var(psi2))/2
		obs[i,] = c(sum(I1, I2, S1, S2), sigma2)
	}
	idx = which(obs[,2]>0)
	log.rc = log(obs[idx,1])
	log.sgm = log(obs[idx,2])
	fit = stats::loess(log.sgm~log.rc)
	df = data.frame(log.rc=log.rc, log.sgm=log.sgm)
	p = ggplot2::ggplot(df, ggplot2::aes(x=log.rc, y=log.sgm)) + ggplot2::geom_point() #+ ggplot2::stat_smooth(method='loess')
	return(list(fit=fit, p=p))
}

estim_group_var.posterior = function(this_data, fit, prior.alpha=1,
	inc_eff_len=2, skp_eff_len=1)
# for inv.gamma distribution, variance is estiamted from 2*alpha
# samples, with sum of squares = 2*beta,
# mode = beta / (alpha+1)
{
	cond_idx = which(this_data$j==1)
	psi1 = this_data$Inc[cond_idx]/inc_eff_len / (this_data$Inc[cond_idx]/inc_eff_len + this_data$Skp[cond_idx]/skp_eff_len)
	psi2 = this_data$Inc[-cond_idx]/inc_eff_len / (this_data$Inc[-cond_idx]/inc_eff_len + this_data$Skp[-cond_idx]/skp_eff_len)
	
	rc = sum(this_data$Inc) + sum(this_data$Skp)
	obs.sgm = (var(psi1) + var(psi2))/2
	prior.sgm = exp(predict(fit, log(rc)))
	obs.alpha = nrow(this_data)/2
	obs.beta = obs.alpha * obs.sgm
	prior.beta = prior.alpha * prior.sgm
	
	post.sgm = (obs.beta+prior.beta ) / (obs.alpha + prior.alpha)
	return(sqrt(post.sgm))
}