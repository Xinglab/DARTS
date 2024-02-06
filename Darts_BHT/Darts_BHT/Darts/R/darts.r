# unreplicated or replicated Bayes factor and Gibbs sampling in Darts setup
# Zijun Zhang, 4.7.2017
# revised 4.27.2017: add paired mode
# revised 5.1.2017: changed the variable name "sigma" in `myLaplace` that could potentially cause a bug
# revised 5.6.2017: add simulator for no replicates
# revised 5.10.2017: changed to true effective length
# revised 5.18.2017: check the sign of determinant in laplace
# revised 7.14.2017: building into a package
# revised 7.17.2017: modified the init value in replicate model to avoid out of domain caused optimization failure
# revised 7.30.2017: modified `test_sim_rdarts` to return ggplot2 AUROC and AUPR curves
# revised 7.30.2017: moved testing function out to a new file `darts_test.r`
# revised 12.1.2017: modified proposal_width in `rdarts_bayesian_gibbs` to adjust for large counts (>10^4)
# revised 12.10.2017: modified init_value and optim_method for better avoiding optim failure in replicate model

#library(pROC)
#library(Rcpp)
#library(MASS)
#library(PRROC)

#sourceCpp('bayes_darts_sampler.cpp')
#sourceCpp('rdarts_sampler.cpp')

## deprecated; using actually effective lengths
# global variables
#inc_eff_len = 2
#skp_eff_len = 1


darts_likelihood = function(mu, delta, data, loglik=F, inc_eff_len=2, skp_eff_len=1)
# likelihood function of darts
# faster version in Rcpp
{
	I1=data$I1; S1=data$S1; I2=data$I2; S2=data$S2
	psi1=mu
	psi2=mu+delta
	invalid = ifelse(loglik, -10^8, 0)
	if(psi2>1 || psi2<0) return(invalid)
	if(psi1>1 || psi1<0) return(invalid)
	eps=0.01
	psi1=ifelse(psi1>1-eps, 1-eps, psi1)
	psi1=ifelse(psi1<eps, eps, psi1)
	psi2=ifelse(psi2>1-eps, 1-eps, psi2)
	psi2=ifelse(psi2<eps, eps, psi2)
	p1 = inc_eff_len*psi1/(inc_eff_len*psi1+skp_eff_len*(1-psi1))
	p2 = inc_eff_len*psi2/(inc_eff_len*psi2+skp_eff_len*(1-psi2))
	if(loglik) {
		return( I1*log(p1) + S1*log(1-p1) + I2*log(p2) + S2*log(1-p2) )
	}
	else{
		return(p1^I1 * (1-p1)^S1 * p2^I2 * (1-p2)^S2)
	}
}



darts_post_pr = function(par, data, tau, loglik=F, inc_eff_len=2, skp_eff_len=1)
# The posterior density of MATS -- for Laplace approximation
# Note this version has no normalizing constant; the normalizing constant
# is to be found by Laplace approximation.
# Prior is a mean-0 variance-tau normal, as a regularization.
{
	mu = par[1]
	delta= par[2]
	#inc_eff_len = 2
	#skp_eff_len = 1
	I1=data$I1; S1=data$S1; I2=data$I2; S2=data$S2
	psi1=mu
	psi2=mu+delta
	invalid = ifelse(loglik, -10^8, 0)
	eps = 10^-2
	if(psi2>1 || psi2<0) return(invalid)
	if(psi1>1 || psi1<0) return(invalid)
	psi1=ifelse(psi1>1-eps, 1-eps, psi1)
	psi1=ifelse(psi1<eps, eps, psi1)
	psi2=ifelse(psi2>1-eps, 1-eps, psi2)
	psi2=ifelse(psi2<eps, eps, psi2)
	p1 = inc_eff_len*psi1/(inc_eff_len*psi1+skp_eff_len*(1-psi1))
	p2 = inc_eff_len*psi2/(inc_eff_len*psi2+skp_eff_len*(1-psi2))
	prior_norm_const = pnorm(1, 0, tau) - pnorm(-1,0,tau)
	if(loglik) {
		lik = I1*log(p1) + S1*log(1-p1) + I2*log(p2) + S2*log(1-p2)
		prior = -delta^2 / (2*tau^2) -log(sqrt(2*pi)*tau) - log(prior_norm_const)
		return(lik + prior)
	}
	else{
		lik = p1^I1 * (1-p1)^S1 * p2^I2 * (1-p2)^S2
		prior = 1/(sqrt(2*pi*tau^2))*exp(-delta^2 / (2*tau^2)) / prior_norm_const
		return(lik * prior)
	}
}




rdarts_post_pr = function(par, data, tau=0.29, sigma=0.05, is_paired=FALSE, inc_eff_len=2, skp_eff_len=1, verbose=FALSE, ...)
# replicate posterior probability
# currently only in log space
{
	n_par = length(par)
	if(is_paired) {
		n_sam = length(unique(data$k))
	} else {
		n_sam = nrow(data)
	}
	mu=par[1]
	delta=par[2]
	mu_k = par[3:(3+n_sam-1)]
	
	invalid = -10^12
	if(any(mu_k>1) || any(mu_k<0) || 
		any((mu_k+delta)>1) || any((mu_k+delta)<0) ||
		mu>1 || mu<0 || mu+delta>1 || mu+delta<0)
		return(invalid)
	
	lik_sum = 0
	# prior likelihood
	prior_norm_const = pnorm(1, 0, tau) - pnorm(-1,0,tau)
	prior = -delta^2 / (2*tau^2) -log(sqrt(2*pi)*tau) - log(prior_norm_const)
	lik_sum = lik_sum + prior
	if(verbose) cat('prior-lik =', lik_sum, '\n')
	
	# group-level likelihood
	lik_sum = lik_sum + sum(dnorm(mu_k, mu, sigma, log=T))
	if(verbose) cat('group-lik =', lik_sum, '\n')
	
	# individual-level likelihood
	for(i in 1:nrow(data))
	{
		Inc = data$Inc[i]
		Skp = data$Skp[i]
		#psi_k = mu_k[i]
		psi_k = mu_k[data$k[i]]
		if(data$j[i]==2)
		{
			psi_k = psi_k + delta
		}
		psi_k = ifelse(psi_k>0.999, 0.999, psi_k)
		psi_k = ifelse(psi_k<0.001, 0.001, psi_k)
		p_k = inc_eff_len*psi_k/(inc_eff_len*psi_k+skp_eff_len*(1-psi_k))
		#lik_sum = lik_sum + dbinom(Inc, Inc+Skp, p_k, log=T)
		lik_sum = lik_sum + Inc*log(p_k) + Skp*log(1-p_k)
		if(verbose) print(Inc*log(p_k) + Skp*log(1-p_k))
	}
	if(verbose) cat('likelihood =', lik_sum, '\n')
	return(lik_sum)
}


rdarts_post_pr_der = function(par, data, tau, sigma, is_paired=FALSE, inc_eff_len=2, skp_eff_len=1, epsilon=1)
{
	#n_sam = nrow(data)
	if(is_paired) {
		n_sam = length(unique(data$k))
	} else {
		n_sam = nrow(data)
	}
	mu = par[1]
	delta = par[2]
	mu_k = par[3:(3+n_sam-1)]
	
	j2_idx = which(data[,'j']==2)
	j1_idx = which(data[,'j']==1)
	if(is_paired)
	{
		stopifnot(length(j1_idx)==length(j2_idx))
		psi_1k = mu_k
		psi_2k = mu_k + delta
	} else {
		psi_1k = mu_k[j1_idx]
		psi_2k = mu_k[j2_idx] + delta
	}
	I_1k = data[j1_idx, 'Inc']
	I_2k = data[j2_idx, 'Inc']
	S_1k = data[j1_idx, 'Skp']
	S_2k = data[j2_idx, 'Skp']
	psi_all = rep(NA, nrow(data))
	psi_all[j1_idx] = psi_1k
	psi_all[j2_idx] = psi_2k
	
	a = inc_eff_len
	b = skp_eff_len
	p_all = psi_all*a / (a*psi_all + b*(1-psi_all))
	p_1 = p_all[j1_idx]
	p_2 = p_all[j2_idx]
	dp_dpsi = a*b/(b*(psi_all-1)-a*psi_all)^2
	
	
	mu_der = sum(-(mu-mu_k))/sigma^2
	#delta_der = -1/tau^2*delta + sum( I_2k/(psi_2k*(1+psi_2k)) - S_2k*2/((1-psi_2k)*(1+psi_2k)) )
	delta_der = -1/tau^2*delta + sum( I_2k/p_2 * dp_dpsi[j2_idx] - S_2k/(1-p_2)*dp_dpsi[j2_idx] )
	

	individual_mu_der = rep(NA, n_sam)
	for(k in 1:n_sam)
	{
		this_idx = which(data$k==k)
		#individual_mu_der[k] = sum(
		#	data[this_idx,'Inc']/(psi_all[this_idx]*(1+psi_all[this_idx])) - 
		#	data[this_idx,'Skp']*2/((1-psi_all[this_idx])*(1+psi_all[this_idx]))
		#	)
		
		individual_mu_der[k] = sum(
			data[this_idx,'Inc']/p_all[this_idx]*dp_dpsi[this_idx] - 
			data[this_idx,'Skp']/(1-p_all[this_idx])*dp_dpsi[this_idx]	
			)
		
	}
	mu_k_der = 1/sigma^2*(mu-mu_k) + individual_mu_der

	
	der = c(mu_der, delta_der, mu_k_der)
	der = der*epsilon
	#print(par)
	#print(der)
	return(der)
}


myLaplace = function(log_func, gr_func, init_value, data, method='Nelder-Mead', logscale=F, verbose=T, ...)
# A very generic implementation of Laplace approximation for integral of posterior density
# log_func: the function pointer to the log posterior density
# data: the data to be passed to log_func
# n_par: number of parameters in log_func
# logscale: return the final integral in log-scale; especially useful for 
#           computing bayes factor and other ratio of two models.
# ...: other keyword arguments passed to log_func
{
	n_par = length(init_value)
	msg = 'None'
	if(method=='L-BFGS-B'){
		res=try(optim(
			par=init_value, 
			fn=log_func,
			gr=gr_func,
			data=data, 
			epsilon=0.001,
			..., 
			hessian=T, 
			control=list(fnscale=-1),
			method='L-BFGS-B',
			lower=c(0.001, -0.999, rep(0.001,length(init_value)-2)),
			upper=c(0.999, 0.999, rep(0.999,length(init_value)-2))
			), silent=T)
		if(inherits(res, "try-error") || sum( abs(res$par - init_value ))<0.001 ) {
			method='Nel'
		} else {
			res$hessian = estim_hess(res$par, gr_func,data, ...)
		}
	}
	if(method!="L-BFGS-B")
	{
		res=optim(
			par=init_value, 
			fn=log_func,
			gr=gr_func,
			data=data, 
			..., 
			hessian=T, 
			control=list(fnscale=-1),
			method=method
			)
	} 
	
	if(verbose) print(res)
	theta_hat = res$par
	hess = res$hessian
	hess[which(is.na(hess),arr.ind=T)]=0
	if(n_par>1)	lap_sigma = try(solve(-hess), silent=T)
	if(n_par==1) lap_sigma = -1/res$hessian
	if(inherits(lap_sigma, "try-error"))
	{
		return(list(MAP=theta_hat, hess=hess, integral=NA))
		#sigma = ginv(-hess)
	}
	det_lap_sigma = ifelse(n_par>1, det(lap_sigma), lap_sigma)
	if(det_lap_sigma<0) {det_lap_sigma=abs(det_lap_sigma);msg='Sigma not PD.'}
	sqrt_sigma = sqrt(det_lap_sigma)
	#if(sqrt_sigma<10^-16 || is.na(sqrt_sigma)) {sqrt_sigma=10^-16; msg='Truncated sqrt_sigma'}
	if(logscale)
	{
		integral=log_func(theta_hat, data, ...) + log(sqrt_sigma) + n_par/2*log(2*pi)
	} else{
		integral=exp(log_func(theta_hat, data, ...)) * sqrt_sigma * (2*pi)^(n_par/2)
	}
	return(list(MAP=theta_hat, hess=hess, integral=integral, msg=msg))
}


rdarts_laplace_bayes_factor = function(data, tau1, tau2, inc_eff_len=2, skp_eff_len=1, sigma=NA, is_paired=F, verbose=FALSE, optim_method="Nel")
# laplace approximation of bayes factor for replicate model
{
	cond_idx = which(data$j==1)
	psi1 = data$Inc[cond_idx]/inc_eff_len / (data$Inc[cond_idx]/inc_eff_len + data$Skp[cond_idx]/skp_eff_len)
	psi1_collapse = sum(data$Inc[cond_idx])/inc_eff_len / sum(data$Inc[cond_idx]/inc_eff_len + data$Skp[cond_idx]/skp_eff_len)
	psi2 = data$Inc[-cond_idx]/inc_eff_len / (data$Inc[-cond_idx]/inc_eff_len + data$Skp[-cond_idx]/skp_eff_len)
	psi2_collapse = sum(data$Inc[-cond_idx])/inc_eff_len / sum(data$Inc[-cond_idx]/inc_eff_len + data$Skp[-cond_idx]/skp_eff_len)
	#mu_mle = mean(psi1)
	mu_mle = psi1_collapse
	mu_mle = ifelse(mu_mle>0.99, 0.99, mu_mle)
	mu_mle = ifelse(mu_mle<0.01, 0.01, mu_mle)
	#delta_mle = mean(psi2)-mean(psi1)
	delta_mle = psi2_collapse - psi1_collapse
	delta_mle = ifelse(mu_mle+delta_mle>0.99, 0.99-mu_mle, delta_mle)
	delta_mle = ifelse(mu_mle+delta_mle<0.01, mu_mle-0.01, delta_mle)
	if(is_paired)
	{
		mu_k_mle = rep(mu_mle, length(unique(data$k)))
	} else {
		mu_k_mle = rep(mu_mle, nrow(data))
	}
	if(is.na(sigma))
		sigma = (sd(psi1) + sd(psi2))/2
	sigma = ifelse(sigma<0.01, 0.01, sigma)
	init_value = c(mu_mle, delta_mle, mu_k_mle)

	m1=myLaplace(rdarts_post_pr, rdarts_post_pr_der, init_value=init_value, data=data, 
		tau=tau1, sigma=sigma, is_paired=is_paired, 
		inc_eff_len=inc_eff_len, skp_eff_len=skp_eff_len,
		logscale=T, verbose=verbose, method=optim_method)

	m2=myLaplace(rdarts_post_pr, rdarts_post_pr_der, init_value=init_value, data=data, 
		tau=tau2, sigma=sigma, is_paired=is_paired, 
		inc_eff_len=inc_eff_len, skp_eff_len=skp_eff_len,
		logscale=T, verbose=verbose, method=optim_method)	
	#if(is.na(m1$integral)||is.na(m2$integral))
	#{
	#	bf=trunc_dnorm_cpp(delta_mle,0,tau1) / trunc_dnorm_cpp(delta_mle,0,tau2)
	#}
	#else
	#{
		bf = exp(m1$integral - m2$integral)
	#}
	return(list(bf=bf, M1_MAP=m1$MAP, M2_MAP=m2$MAP))
}



darts_laplace_bayes_factor = function(tau1, tau2, data, 
	inc_eff_len=2, skp_eff_len=1, verbose=F, optim_method='Nel')
# the Laplace approximation version of bf
{
	I1 = data$I1
	S1 = data$S1
	I2 = data$I2
	S2 = data$S2
	mu_init = I1/inc_eff_len /(I1/inc_eff_len + S1/skp_eff_len)
	delta_init = I2/inc_eff_len /(I2/inc_eff_len + S2/skp_eff_len) - mu_init
	
	mu_init = ifelse(mu_init>0.99, 0.99, mu_init)
	mu_init = ifelse(mu_init<0.01, 0.01, mu_init)
	delta_init = ifelse(mu_init+delta_init>0.99, 0.99-mu_init, delta_init)
	delta_init = ifelse(mu_init+delta_init<0.01, mu_init-0.01, delta_init)
	
	init_value = c(mu_init, delta_init)
	M1=myLaplace(darts_post_pr, darts_post_der, init_value, data, 
		inc_eff_len=inc_eff_len, skp_eff_len=skp_eff_len,
		logscale=T, method=optim_method, verbose=verbose, tau=tau1, loglik=T)
	M2=myLaplace(darts_post_pr, darts_post_der, init_value, data, 
		inc_eff_len=inc_eff_len, skp_eff_len=skp_eff_len,
		logscale=T, method=optim_method, verbose=verbose, tau=tau2, loglik=T)	
	if(is.na(M1$integral)||is.na(M2$integral))
	{
		bf=trunc_dnorm_cpp(delta_init,0,tau1) / trunc_dnorm_cpp(delta_init,0,tau2)
	} else {
		bf=exp(M1$integral - M2$integral)
	}
	return(list(bf=bf, M1_MAP=M1$MAP, M2_MAP=M2$MAP))
}


rdarts_bayesian_inference = function(rho, data, tau1=0.28, tau2=0.02,
	inc_eff_len=2, skp_eff_len=1, sigma=NA, C=0.05, N=5000, is_paired=F, optim_method='Nel')
{
	# handling missing sigma
	cond_idx = which(data$j==1)
	cond_idx = which(data$j==1)
	psi1 = data$Inc[cond_idx]/inc_eff_len / (data$Inc[cond_idx]/inc_eff_len + data$Skp[cond_idx]/skp_eff_len)
	psi2 = data$Inc[-cond_idx]/inc_eff_len / (data$Inc[-cond_idx]/inc_eff_len + data$Skp[-cond_idx]/skp_eff_len)
	if(is.na(sigma))
		sigma = (sd(psi1) + sd(psi2))/2
	
	# now actual inference..	
	bf_res = rdarts_laplace_bayes_factor(data, tau1, tau2,
		inc_eff_len=inc_eff_len, skp_eff_len=skp_eff_len,
		sigma=sigma, is_paired=is_paired, optim_method=optim_method)
	bf = bf_res$bf
	m1_init_val = bf_res$M1_MAP
	m2_init_val = bf_res$M2_MAP
	## NOTE: this can cause problems in the future:
	bf = ifelse(is.na(bf), 1, bf)
	bf = ifelse(bf>1000, 1000, bf)  # need to bound bf to avoid wild behaviours
	bf = ifelse(bf<0.001, 0.001, bf)
	a = rho / (1-rho) * bf
	pr = a / (a+1)
	
	# for large counts, need smaller width to sample/explore the highly-concentrated probability space
	# added 12.1.2017
	proposal_width = ifelse(sum(data)>10^4, 0.01, 0.02) 
	m1_post = rdarts_posterior_MCMC_sampler(N*10+100, tau1, sigma, as.matrix(data), m1_init_val, inc_eff_len=inc_eff_len, skp_eff_len=skp_eff_len, burnin=100, thinning=10, is_paired=is_paired, proposal_width=proposal_width)
	m2_post = rdarts_posterior_MCMC_sampler(N*10+100, tau2, sigma, as.matrix(data), m2_init_val, inc_eff_len=inc_eff_len, skp_eff_len=skp_eff_len, burnin=100, thinning=10, is_paired=is_paired, proposal_width=proposal_width)
	
	post = matrix(NA, ncol=length(m1_init_val), nrow=N)
	if(is.na(pr)) return(list(post_pr=NA, mix_pr=pr, post=post))
	for(i in 1:N)
	{
		if(runif(1)<pr) { # H1/M1 selected
			post[i,] = m1_post[sample(1:N, 1),]
		} else { # H0/M2 selected
			post[i,] = m2_post[sample(1:N,1),]
		}
	}
	post_pr = mean(abs(post[,2])>C)
	return(list(post_pr=post_pr, mix_pr=pr, post=post))
}


estim_hess = function(theta_hat, gr_func, data, eps=0.0005, ...)
{
	n_par = length(theta_hat)
	hess = matrix(0, nrow=n_par, ncol=n_par)
	for(i in 1:n_par)
	{
		#gr_i = gr_func(theta_hat, data, epsilon=1, ...)[i]
		for(j in 1:n_par)
		{
			new_theta_1 = theta_hat
			new_theta_2 = theta_hat
			new_theta_1[j] = new_theta_1[j]+eps
			new_theta_2[j] = new_theta_2[j]-eps
			gr_i_dj_1 = gr_func(new_theta_1, data, epsilon=1, ...)[i]
			gr_i_dj_2 = gr_func(new_theta_2, data, epsilon=1, ...)[i]
			hess_i_j = (gr_i_dj_1 - gr_i_dj_2)/(2*eps)
			hess[i,j] = hess_i_j
			hess[j,i] = hess_i_j
		}
	}
	return(hess)
}


rdarts_bayesian_gibbs = function(rho, data, tau1=0.28, tau2=0.02, sigma=NA, 
	N=1500, burnin=100, thinning=10, is_paired=F, sgm_prior_fit=NA,
	inc_eff_len=2, skp_eff_len=1, optim_method='L-BFGS-B')
# for calling from rpy2
{
	# handling missing sigma
	if(class(sgm_prior_fit)!="loess") {
		cond_idx = which(data$j==1)
		psi1 = data$Inc[cond_idx]/inc_eff_len / (data$Inc[cond_idx]/inc_eff_len + data$Skp[cond_idx]/skp_eff_len)
		psi2 = data$Inc[-cond_idx]/inc_eff_len / (data$Inc[-cond_idx]/inc_eff_len + data$Skp[-cond_idx]/skp_eff_len)
		if(is.na(sigma))
			sigma = (sd(psi1) + sd(psi2))/2
		sigma = ifelse(sigma<0.01, 0.01, sigma)
	} else {
		prior.alpha = 2
		sigma = estim_group_var.posterior(data, sgm_prior_fit, prior.alpha=prior.alpha, 
			inc_eff_len=inc_eff_len, skp_eff_len=skp_eff_len)
	}
	# now actual inference..	
	bf_res = rdarts_laplace_bayes_factor(data, tau1, tau2,
		inc_eff_len=inc_eff_len, skp_eff_len=skp_eff_len,
		sigma=sigma, is_paired=is_paired, optim_method=optim_method)
	bf = bf_res$bf
	m1_init_val = bf_res$M1_MAP
	m2_init_val = bf_res$M2_MAP
	
	## NOTE: this can cause problems in the future:
	bf = ifelse(is.na(bf), 1, bf)
	bf = ifelse(bf>1000, 1000, bf)  # need to bound bf to avoid wild behaviours
	bf = ifelse(bf<0.001, 0.001, bf)
	a = rho / (1-rho) * bf
	pr = a / (a+1)
	
	proposal_width = ifelse(sum(data)>10^4, 0.01, 0.02) 
	m1_post = rdarts_posterior_MCMC_sampler(N*thinning+burnin, tau1, sigma, cbind(data$k, data$j, data$Inc, data$Skp), m1_init_val, burnin=burnin, thinning=thinning, is_paired=is_paired, proposal_width=proposal_width)
	m2_post = rdarts_posterior_MCMC_sampler(N*thinning+burnin, tau2, sigma, cbind(data$k, data$j, data$Inc, data$Skp), m2_init_val, burnin=burnin, thinning=thinning, is_paired=is_paired, proposal_width=proposal_width)
	
	post = matrix(NA, ncol=length(m1_init_val), nrow=N)
	if(is.na(pr)) return(list(post_pr=NA, mix_pr=pr, post=post))
	for(i in 1:N)
	{
		if(runif(1)<pr) { # H1/M1 selected
			post[i,] = m1_post[sample(1:N, 1),]
		} else { # H0/M2 selected
			post[i,] = m2_post[sample(1:N,1),]
		}
	}
	return(post)
}


darts_bayesian_gibbs = function(rho, data, tau1=0.5, tau2=0.045, C=0.05, verbose=FALSE,
	N=1500, burnin=100, thinning=10, inc_eff_len=2, skp_eff_len=1)
# # for calling from rpy2
{
	bf_res = darts_laplace_bayes_factor(tau1, tau2, data, 
		inc_eff_len=inc_eff_len, skp_eff_len=skp_eff_len, optim_method='Nel')
	bf = bf_res$bf
	m1_init_val = bf_res$M1_MAP
	m2_init_val = bf_res$M2_MAP
	bf = ifelse(bf>1000, 1000, bf)  # need to bound bf to avoid wild behaviours
	a = rho / (1-rho) * bf
	pr = a / (a+1)
	
	m1_post = darts_posterior_MCMC_sampler(N*thinning+burnin, tau1, data, m1_init_val, burnin=burnin, thinning=thinning, inc_eff_len=inc_eff_len, skp_eff_len=skp_eff_len)
	m2_post = darts_posterior_MCMC_sampler(N*thinning+burnin, tau2, data, m2_init_val, burnin=burnin, thinning=thinning, inc_eff_len=inc_eff_len, skp_eff_len=skp_eff_len)

	post = matrix(NA, ncol=2, nrow=N)
	#if(verbose) cat('sampling joint-model posterior..\n')
	if(is.na(pr)) return(post)
	for(i in 1:N)
	{
		if(runif(1)<pr) { # H1/M1 selected
			post[i,] = m1_post[sample(1:N, 1),]
		} else { # H0/M2 selected
			post[i,] = m2_post[sample(1:N,1),]
		}
	}
	return(post)
}


darts_bayesian_inference = function(rho, data, tau1=0.5, tau2=0.045, C=0.05, N=5000, verbose=FALSE,
	inc_eff_len=2, skp_eff_len=1)
# draw posterior samples to do inference
{
	bf_res = darts_laplace_bayes_factor(tau1, tau2, data, 
		inc_eff_len=inc_eff_len, skp_eff_len=skp_eff_len, optim_method='Nel')
	bf = bf_res$bf
	m1_init_val = bf_res$M1_MAP
	m2_init_val = bf_res$M2_MAP
	bf = ifelse(bf>100, 100, bf)  # need to bound bf to avoid wild behaviours
	bf = ifelse(bf<0.01, 0.01, bf)  # need to bound bf to avoid wild behaviours
	a = rho / (1-rho) * bf
	pr = a / (a+1)
	#if(verbose) cat('sampling model-wise posterior..\n')
	
	m1_post = darts_posterior_MCMC_sampler(N*10+100, tau1, data, m1_init_val, burnin=100, thinning=10, inc_eff_len=inc_eff_len, skp_eff_len=skp_eff_len)
	m2_post = darts_posterior_MCMC_sampler(N*10+100, tau2, data, m2_init_val, burnin=100, thinning=10, inc_eff_len=inc_eff_len, skp_eff_len=skp_eff_len)
	#####
	# slower version of rejection sampling
	# deprecated, but can serve as baseline for 
	# checking consistency
	#m1_post = darts_rej_samp(N, tau1, data)
	#m2_post = darts_rej_samp(N, tau2, data)
	#####
	
	post = matrix(NA, ncol=2, nrow=N)
	#if(verbose) cat('sampling joint-model posterior..\n')
	if(is.na(pr)) return(post)
	for(i in 1:N)
	{
		if(runif(1)<pr) { # H1/M1 selected
			post[i,] = m1_post[sample(1:N, 1),]
		} else { # H0/M2 selected
			post[i,] = m2_post[sample(1:N,1),]
		}
	}
	post_pr = mean(abs(post[,2])>C)
	return(list(post_pr=post_pr, mix_pr=pr, post=post))
}


darts_post_der = function(par, data, tau, inc_eff_len=2, skp_eff_len=1, loglik=T, epsilon=1)
# the first derivative of darts posterior density
# currently only implemented loglik form
# seems still problemtic... DEPRECATED as of 2.7.2017
# re-used as of 5.11.2017
{
	mu = par[1]
	delta = par[2]
	#inc_eff_len = 2
	#skp_eff_len = 1
	I1 = data$I1
	S1 = data$S1
	I2 = data$I2
	S2 = data$S2
	psi1 = mu; psi2 = mu+delta
	p1 = inc_eff_len*psi1/(inc_eff_len*psi1+skp_eff_len*(1-psi1))
	p2 = inc_eff_len*psi2/(inc_eff_len*psi2+skp_eff_len*(1-psi2))
	
	a = inc_eff_len
	b = skp_eff_len
	dp1_dpsi1 = a*b/(b*(psi1-1)-a*psi1)^2
	dp2_dpsi2 = a*b/(b*(psi2-1)-a*psi2)^2
	
	#mu_gr = I1 / (mu*(1+mu)) - S1*2/((1-mu)*(1+mu)) +
	#	    I2 / ((mu+delta)*(1+mu+delta)) - S2*2/((1-mu-delta)*(1+mu+delta))
	mu_gr = I1/p1 * dp1_dpsi1 - S1/(1-p1) * dp1_dpsi1 + 
			I2/p2 * dp2_dpsi2 - S2/(1-p2) * dp2_dpsi2
	
	#delta_gr = I2 / ((mu+delta)*(1+mu+delta)) -
	#           S2*2/((1-mu-delta)*(1+mu+delta)) -
	#		   delta/tau^2
	delta_gr = I2/p2 * dp2_dpsi2 - S2/(1-p2) * dp2_dpsi2 - 
			   delta/tau^2
	
	der = c(mu_gr, delta_gr)
	der = der*epsilon
	return(der)
}



#####
## Command-line usage


