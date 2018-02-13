## parse rMATS to the same set of simulated data
## for comparison
## Zijun Zhang
## 2.3.2018

library(PRROC)

#RMATS_BIN = '/u/nobackup/yxing/NOBACKUP/frankwoe/src/rMATS-STAT/rMATS_unpaired.py'
RMATS_BIN = '..\\rMATS-STAT\\rMATS_unpaired.py'

plot_pr=function(pr_df, eval_df, config_list)
{
	covg = config_list$covg
	K = config_list$K
	n = config_list$n
	outlier = config_list$outlier
	sigma = config_list$sigma
	pr_df = pr_df[pr_df$group!='rDarts_Nel',]
	p_pr = ggplot2::ggplot(pr_df, ggplot2::aes(x=RE, y=PR, group=group, colour=group)) +
		ggplot2::geom_line(size=1) +
		ggplot2::scale_colour_manual(
			name=NULL,
			values=c('black', 'red', 'darkblue'),
			labels=c(paste(c('Pool Darts, aupr=', 'rDarts, aupr=', 'rMATS, aupr='), 
				round(c(eval_df$aupr[2], eval_df$aupr[3], eval_df$aupr[4]),2), sep=''))) +
		ggplot2::theme_bw() + 
		ggplot2::ggtitle(
			paste0(
				'Sim.AUPR: covg=',covg,
				', K=',K,
				', n=',n, ',\n', 
				'outlier=',outlier,
				#', is_paired=',is_paired,
				', sigma=',sigma)) +
			Darts_theme
	return(p_pr)
}


plot_roc=function(roc_df, eval_df, config_list)
{
	covg = config_list$covg
	K = config_list$K
	n = config_list$n
	outlier = config_list$outlier
	sigma = config_list$sigma
	roc_df = roc_df[roc_df$group!='rDarts_Nel',]
	p_roc = ggplot2::ggplot(roc_df, ggplot2::aes(x=FPR, y=TPR, group=group, colour=group)) +
		ggplot2::geom_line(size=1) +
		ggplot2::scale_colour_manual(
			name=NULL,
			values=c('black', 'red', 'darkblue'),
			labels=c(paste(c('Pool Darts, auroc=', 'rDarts, auroc=', 'rMATS, auroc='), 
				round(c(eval_df$auroc[2], eval_df$auroc[3], eval_df$auroc[4]),2), sep=''))) +		
		ggplot2::ggtitle(
			paste0(
				'Sim.AUROC: covg=',covg, 
				', K=',K,
				', n=',n, ',\n', 
				'outlier=',outlier,
				#', is_paired=',is_paired,
				', sigma=',sigma)) +
		Darts_theme
	return(p_roc)
}



parse_rMATS = function(res, config_list)
{
	sim_df = res$sim_df
	write.table(sim_df, file='inc.txt', sep='\t', row.names=F, col.names=T, quote=F)
	#cmd = paste0('python ', RMATS_BIN, ' inc.txt . 4 0.05 1>/dev/null 2>&1')
	cmd = paste0('python ', RMATS_BIN, ' inc.txt . 1 0.05')
	system(cmd, ignore.stdout=T, ignore.stderr=T, wait=T)
	data = read.table('rMATS_Result_P.txt', header=T)
	roc = roc.curve(-data$PValue[abs(data$sim.delta)>0.05], -data$PValue[abs(data$sim.delta)<=0.05], curve=T)
	pr = pr.curve(-data$PValue[abs(data$sim.delta)>0.05], -data$PValue[abs(data$sim.delta)<=0.05], curve=T)
	roc_df = rbind.data.frame(
		res$roc_df,
		data.frame(FPR=roc$curve[,1], TPR=roc$curve[,2], group='rMATS', stringsAsFactors=F),
		stringsAsFactors=F
	)
	pr_df = rbind.data.frame(
		res$pr_df,
		data.frame(PR=pr$curve[,2], RE=pr$curve[,1], group='rMATS', stringsAsFactors=F),
		stringsAsFactors=F
	)
	eval_df = rbind.data.frame(
		res$eval_df,
		data.frame(method='rMATS', auroc=roc$auc, aupr=pr$auc.integral, stringsAsFactors=F), 
		stringsAsFactors=F
	)
	
	new.p_pr = plot_pr(pr_df, eval_df, config_list)
	new.p_roc = plot_roc(roc_df, eval_df, config_list)
	return(list(eval_df=eval_df, p_pr=new.p_pr, p_roc=new.p_roc))
}
