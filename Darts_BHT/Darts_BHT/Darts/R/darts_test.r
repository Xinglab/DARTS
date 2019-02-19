# For testing no-replicate, replicate Darts model using simulation data
# Zijun Zhang
# 7.30.2017
# revised 2.3.2018: output simulation data from `test_sim_rdarts`

Darts_theme = ggplot2::theme(
		legend.text = ggplot2::element_text(size = 13),
		plot.title = ggplot2::element_text(size=13, face="bold"), 
		axis.title.y = ggplot2::element_text(size=13), 
		axis.title.x = ggplot2::element_text(size=13),
		axis.text.y = ggplot2::element_text(size=13, angle = 90, hjust = 0.5, vjust=0.5), 
		axis.text.x = ggplot2::element_text(size=13, angle=0, hjust=0.5, vjust=0.5),
		legend.background = ggplot2::element_rect(fill = "transparent", colour = "transparent")) +
	ggplot2::theme(
		panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
		panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
		ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
		panel.grid.minor = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black"))						


sim_data = function(lambda=0.1, covg=10, C=0.01, inc_eff_len=2, skp_eff_len=1)
# a simple simulator to generate mu, delta and data(i.e. I1,S1,I2,S2)
{
	#inc_eff_len = 2
	#skp_eff_len = 1
	pass=FALSE
	while(!pass)
	{
		mu = runif(1,0.1,0.9)
		if(runif(1)<lambda)
		{ # is a significant change
			delta = runif(1, C, 0.4) * ifelse(runif(1)>0.5, 1, -1)
		} else { # not a significant change
			delta = runif(1, -C, C)
		}
		if(mu+delta>0.01 && mu+delta<0.99) pass=TRUE
	}
	pass.2=FALSE
	while(!pass.2)
	{
		psi1=mu; psi2=mu+delta
		p1 = inc_eff_len*psi1/(inc_eff_len*psi1+skp_eff_len*(1-psi1))
		p2 = inc_eff_len*psi2/(inc_eff_len*psi2+skp_eff_len*(1-psi2))
		#n1 = rpois(1,covg)
		n1 = covg
		I1 = rbinom(1, p=p1, size=n1); S1=n1-I1
		#n2 = rpois(1,covg)
		n2 = covg
		I2 = rbinom(1, p=p2, size=n2); S2=n2-I2
		p1=I1/(I1+S1)
		p2=I2/(I2+S2)
		if(! (p1==p2 && (p1==0 || p1==1) )) pass.2=TRUE
	}
	return(list(mu=mu, delta=delta, inc_eff_len=inc_eff_len, skp_eff_len=skp_eff_len, data=list(I1=I1,S1=S1,I2=I2,S2=S2)))
}


sim_replicate_data = function(K=20, sigma=0.05, covg=50, outlier=FALSE, is_paired=FALSE,
	inc_eff_len=2, skp_eff_len=1)
# simulate read counts for unpaired/paired replicates
{
	if(is_paired)
		K=round(K/2)
	lambda = 0.5 # mixture prob.; higher this value, fewer differential events
	tau1 = 0.28
	tau2 = 0.02
	sim_res = list()
	gen_mdl=0
	while(TRUE)
	{
		mu = runif(1)
		p_lambda = runif(1)
		if(p_lambda < lambda)
		{
			delta = rnorm(1,0,tau2)
			gen_mdl=2
		} else {
			delta = rnorm(1,0,tau1)
			gen_mdl=1
		}
		if(mu+delta<0.99 && mu+delta>0.01 && mu<0.99 && mu>0.01)
			break
	}
	sim_res$mu = mu
	sim_res$delta = delta
	psi1 = mu
	psi2 = mu+delta
	Inc=c(); Skp=c();
	jj=c(); kk=c(); # jj index conditions "1" or "2"; kk index sample
	psi_kk = c()
	mu_kk = c()
	for(k in 1:K)
	{
		while(TRUE)
		{
			if(outlier && k==K) {
				mu_k = mu+rnorm(1,0,sigma*10)
				delta_out = delta + rnorm(1,0,sigma)
				if(mu_k+delta_out<0.99 && mu_k+delta_out>0.01 && mu_k<0.99 && mu_k>0.01)
					break
			} else {
				mu_k = mu + rnorm(1,0,sigma)
				if(mu_k+delta<0.99 && mu_k+delta>0.01 && mu_k<0.99 && mu_k>0.01)
					break
			}
			
		}
		if(outlier && k==K) {
			this_data = sim_individual_data(mu_k, delta_out, covg=covg, inc_eff_len=inc_eff_len, skp_eff_len=skp_eff_len)
		} else {
			this_data = sim_individual_data(mu_k, delta, covg=covg, inc_eff_len=inc_eff_len, skp_eff_len=skp_eff_len)
		}
		if(is_paired)
		{
			Inc = c(Inc, this_data$data$I1)
			Skp = c(Skp, this_data$data$S1)
			psi_k = mu_k
			
			kk = c(kk, k)
			jj = c(jj, 1)
			psi_kk = c(psi_kk, psi_k)
			mu_kk = c(mu_kk, mu_k)
			
			Inc = c(Inc, this_data$data$I2)
			Skp = c(Skp, this_data$data$S2)
			psi_k = mu_k + delta

			kk = c(kk, k)
			jj = c(jj, 2)
			psi_kk = c(psi_kk, psi_k)
			mu_kk = c(mu_kk, mu_k)
		} else {
			if(k%%2)
			{
				Inc = c(Inc, this_data$data$I1)
				Skp = c(Skp, this_data$data$S1)
				psi_k = mu_k
			} else {
				Inc = c(Inc, this_data$data$I2)
				Skp = c(Skp, this_data$data$S2)
				psi_k = mu_k + delta
			}
			kk = c(kk, k)
			jj = c(jj, 2-k%%2)
			psi_kk = c(psi_kk, psi_k)
			mu_kk = c(mu_kk, mu_k)
		}
		
	}
	data = data.frame(k=kk, j=jj, Inc=Inc, Skp=Skp, psi_k=psi_kk, mu_k=mu_kk)
	#data$inc_eff_len=2
	#data$skp_eff_len=1
	return(list(mu=mu, delta=delta, gen_mdl=gen_mdl, 
		inc_eff_len=inc_eff_len, skp_eff_len=skp_eff_len,
		data=data))
}



sim_individual_data = function(mu_k, delta_k, covg=20, inc_eff_len=2, skp_eff_len=1)
# a simple simulator to generate mu, delta and data(i.e. I1,S1,I2,S2)
{
	#inc_eff_len = 2
	#skp_eff_len = 1
	mu = mu_k
	delta = delta_k
	pass.2=FALSE
	while(!pass.2)
	{
		psi1=mu; psi2=mu+delta
		p1 = inc_eff_len*psi1/(inc_eff_len*psi1+skp_eff_len*(1-psi1))
		p2 = inc_eff_len*psi2/(inc_eff_len*psi2+skp_eff_len*(1-psi2))
		n1 = rpois(1,covg)
		#n1 = covg
		I1 = rbinom(1, p=p1, size=n1); S1=n1-I1
		n2 = rpois(1,covg)
		#n2 = covg
		I2 = rbinom(1, p=p2, size=n2); S2=n2-I2
		p1=I1/(I1+S1)
		p2=I2/(I2+S2)
		if(! (p1==p2 && (p1==0 || p1==1) )) pass.2=TRUE
	}
	return(list(mu=mu, delta=delta, data=list(I1=I1,S1=S1,I2=I2,S2=S2)))
}


collapse_replicates = function(data)
{
	I1 = sum(data$Inc[data$j==1])
	I2 = sum(data$Inc[data$j==2])
	S1 = sum(data$Skp[data$j==1])
	S2 = sum(data$Skp[data$j==2])
	return(list(I1=I1, S1=S1, I2=I2, S2=S2))
}



test_sim_rdarts = function(n=300, K=10, sigma=0.05, covg=30, outlier=T, is_paired=FALSE, seed=12345)
{
	set.seed(seed)
	inc_eff_len = 2 
	skp_eff_len = 1
	sim_df = matrix(NA, nrow=n, ncol=9)
	colnames(sim_df) = c('ID', 'I1', 'S1', 'I2', 'S2', 'inc_len', 'skp_len', 'sim.mu', 'sim.delta')
	if(is_paired)
	{
		res = matrix(NA, nrow=n, ncol=5)
		colnames(res) = c('gen_mdl', 'delta', 'post_r', 'post_s', 'post_ur')
	} else {
		#res = matrix(NA, nrow=n, ncol=4)
		#colnames(res) = c('gen_mdl', 'delta', 'post_r', 'post_s')
		res = matrix(NA, nrow=n, ncol=5)
		colnames(res) = c('gen_mdl', 'delta', 'post_r', 'post_s', 'post_bfgs')
	}
	for(i in 1:n)
	{
		if(! i%%10)
			cat('[', format(Sys.time(), "%a %b %d %X %Y"),'] ', i, '\n', sep='')
		sim=sim_replicate_data(K=K, sigma=sigma, covg=covg, outlier=outlier, is_paired=is_paired,
			inc_eff_len=inc_eff_len, skp_eff_len=skp_eff_len)
		data = sim$data
		
		#hist(data$psi_k[data$j==1], ylim=c(0,10), xlim=c(min(data$psi_k)*0.8, max(data$psi_k)*1.2), col=rgb(1,0,0,0.8), freq=T, breaks=15)
		#hist(data$psi_k[data$j==2], xlim=c(min(data$psi_k)*0.8, max(data$psi_k)*1.2), col=rgb(0,0,1,0.8), add=T, freq=T, breaks=15)
		
		post_r = rdarts_bayesian_inference(0.5, data, tau1=0.3, tau2=0.03, sigma=NA, C=0.05, is_paired=is_paired,
			inc_eff_len=inc_eff_len, skp_eff_len=skp_eff_len, optim_method='Nel')
		post_r2 = rdarts_bayesian_inference(0.5, data, tau1=0.3, tau2=0.03, sigma=NA, C=0.05, is_paired=is_paired,
			inc_eff_len=inc_eff_len, skp_eff_len=skp_eff_len, optim_method='L-BFGS-B')
		
		
		if(is_paired)
			post_ur = rdarts_bayesian_inference(0.5, data, tau1=0.3, tau2=0.03, sigma=NA, C=0.05, is_paired=F,
					inc_eff_len=inc_eff_len, skp_eff_len=skp_eff_len)
		
		cdata = collapse_replicates(data)
		post_s = darts_bayesian_inference(0.5, cdata, tau1=0.3, tau2=0.03, C=0.05,
				inc_eff_len=inc_eff_len, skp_eff_len=skp_eff_len)
		
		# 'ID', 'I1', 'S1', 'I2', 'S2', 'inc_len', 'skp_len', 'sim.mu', 'sim.delta'
		sim_df[i,] = c(i, 
			paste0(data[data$j==1,3], collapse=','), 
			paste0(data[data$j==1,4], collapse=','),
			paste0(data[data$j==2,3], collapse=','),
			paste0(data[data$j==2,4], collapse=','),
			inc_eff_len,
			skp_eff_len,
			sim$mu,
			sim$delta)
		if(is_paired) {
			res[i,] = c(sim$gen_mdl, sim$delta, post_r$post_pr, post_s$post_pr, post_ur$post_pr)
		} else {
			res[i,] = c(sim$gen_mdl, sim$delta, post_r$post_pr, post_s$post_pr, post_r2$post_pr)
			#res[i,] = c(sim$gen_mdl, sim$delta, post_r2$post_pr, post_s$post_pr)
		}
	}
	res=as.data.frame(res)
	
	roc_1=PRROC::roc.curve(res$post_r[abs(res$delta)>0.05], res$post_r[abs(res$delta)<=0.05], curve=T)
	roc_2=PRROC::roc.curve(res$post_s[abs(res$delta)>0.05], res$post_s[abs(res$delta)<=0.05], curve=T)
	roc_3=PRROC::roc.curve(res$post_bfgs[abs(res$delta)>0.05], res$post_bfgs[abs(res$delta)<=0.05], curve=T)
	
	roc_df = rbind.data.frame(
		data.frame(FPR=roc_1$curve[,1], TPR=roc_1$curve[,2], group='rDarts_Nel', stringsAsFactors=F),
		data.frame(FPR=roc_2$curve[,1], TPR=roc_2$curve[,2], group='Darts', stringsAsFactors=F),
		data.frame(FPR=roc_3$curve[,1], TPR=roc_3$curve[,2], group='rDarts_L-BFGS-B', stringsAsFactors=F),	
		stringsAsFactors=F)
	p_roc = ggplot2::ggplot(roc_df, ggplot2::aes(x=FPR, y=TPR, group=group, colour=group)) +
		ggplot2::geom_line(size=1) +
		ggplot2::scale_colour_manual(
			name=NULL,
			values=c('black', 'red', 'darkblue'),
			labels=c(paste(c('Pool Darts, auc=', 'rDarts_L-BFGS-B, auc=', 'rDarts_Nel, auc='), round(c(roc_2$auc, roc_3$auc, roc_1$auc),2), sep=''))) +
		ggplot2::theme_bw() + 
		ggplot2::ggtitle(
			paste0(
				'Sim.AUROC: covg=',covg, 
				', K=',K,
				', n=',n, 
				', outlier=',outlier, '\n',
				', is_paired=',is_paired,
				', sigma=',sigma)) +
		Darts_theme
	
	pr_1 = PRROC::pr.curve(res$post_r[abs(res$delta)>0.05], res$post_r[abs(res$delta)<=0.05], curve=T)
	pr_2 = PRROC::pr.curve(res$post_s[abs(res$delta)>0.05], res$post_s[abs(res$delta)<=0.05], curve=T)
	pr_3 = PRROC::pr.curve(res$post_bfgs[abs(res$delta)>0.05], res$post_bfgs[abs(res$delta)<=0.05], curve=T)
	
	pr_df = rbind.data.frame(
		data.frame(PR=pr_1$curve[,2], RE=pr_1$curve[,1], group='rDarts_Nel', stringsAsFactors=F),
		data.frame(PR=pr_2$curve[,2], RE=pr_2$curve[,1], group='Darts', stringsAsFactors=F),
		data.frame(PR=pr_3$curve[,2], RE=pr_3$curve[,1], group='rDarts_L-BFGS-B', stringsAsFactors=F),		
		stringsAsFactors=F)
		
	p_pr = ggplot2::ggplot(pr_df, ggplot2::aes(x=RE, y=PR, group=group, colour=group)) +
		ggplot2::geom_line(size=1) +
		ggplot2::scale_colour_manual(
			name=NULL,
			values=c('black', 'red', 'darkblue'),
			labels=c(paste(c('Pool Darts, aupr=', 'rDarts_L-BFGS-B, aupr=', 'rDarts_Nel, aupr='), round(c(pr_2$auc.integral, pr_3$auc.integral, pr_1$auc.integral),2), sep=''))) +
		ggplot2::theme_bw() + 
		ggplot2::ggtitle(
			paste0(
				'Sim.AUPR: covg=',covg,
				', K=',K,
				', n=',n, 
				', outlier=',outlier, '\n',
				', is_paired=',is_paired,
				', sigma=',sigma)) +
		Darts_theme
	
	eval_df = data.frame(method=c('Nel','Pool','L-BFGS-B'), auroc=c(roc_1$auc, roc_2$auc, roc_3$auc), 
		aupr=c(pr_1$auc.integral, pr_2$auc.integral, pr_3$auc.integral))
		
	print(eval_df)
	sim_df = as.data.frame(sim_df)
	return(list(sim_df=sim_df, eval_df=eval_df, p_roc = p_roc, p_pr = p_pr, roc_df=roc_df, pr_df=pr_df))
	
}

test_sim_darts = function(covg=50)
{
	sim = sim_data(lambda=0.5, covg=covg, inc_eff_len=1, skp_eff_len=1)
	res = darts_bayesian_inference(rho=0.5, data=sim$data, tau1=0.3, tau2=0.03, inc_eff_len=1, skp_eff_len=1)
	return(list(sim=sim, res=res$post_pr))
}
 