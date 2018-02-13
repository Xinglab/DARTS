## plot the comparison of repDarts in a single big pdf
# Zijun Zhang
# 8.2.2017
# revised 2.3.2018: added rMATS parser

library(Darts)
library(ggplot2)
library(gridExtra)

source('sim_rmats_parser.R')

config = data.frame(
	K    = c( 6, 6, 6, 10, 10, 10 ),
	sgm  = c( 0.05, 0.05, 0.35, 0.05, 0.05, 0.35 ),
	covg = c( 50, 50, 50, 30, 30, 30 ),
	outlr= c( F, T, F, F, T, F),
	n    = 3000,
	pair = c( F, F, F, F, F, F ) )


# config = data.frame(
	# K    = c( 6, 6, 6),
	# sgm  = c( 0.05, 0.05, 0.35),
	# covg = c( 50, 50, 50),
	# outlr= c( F, T, F),
	# n    = 1000,
	# pair = c( F, F, F) )
	
	
roc_list = list()
pr_list = list()

for(i in 1:nrow(config))
{
	res = test_sim_rdarts(
		n=config$n[i],
		K=config$K[i],
		sigma=config$sgm[i],
		covg=config$covg[i],
		outlier=config$outlr[i],
		is_paired=config$pair[i],
		seed=12345
		)
	config_list = list(
		n=config$n[i],
		K=config$K[i],
		sigma=config$sgm[i],
		covg=config$covg[i],
		outlier=config$outlr[i]
		)
	new_res = parse_rMATS(res, config_list)
	roc_list[[i]] = new_res$p_roc
	pr_list[[i]] = new_res$p_pr
	plot(new_res$p_pr)
}

pdf('roc_all.2row.pdf', useDingbats=F,
	width=15, height=6)
do.call(grid.arrange, c(roc_list, list(ncol=3)))
dev.off()


pdf('pr_all.2row.pdf', useDingbats=F,
	width=15, height=6)
do.call(grid.arrange, c(pr_list, list(ncol=3)))
dev.off()