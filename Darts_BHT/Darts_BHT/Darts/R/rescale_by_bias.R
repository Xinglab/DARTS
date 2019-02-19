## a new rescaling scheme based on 
## estimation of a new bias term
## Zijun Zhang
## 3.19.2018


rescale_by_bias = function(y_pred, y_true, b_inherent=0.05554411)
{
	nb_pos = sum(y_true>0.9, na.rm=T)
	nb_neg = sum(y_true<0.1, na.rm=T)
	log_ratio = log(y_pred) - log(1-y_pred)
	log_ratio = log_ratio - log(nb_pos/nb_neg) - b_inherent ## inherent bias term
	new_pred_ratio = exp(log_ratio)
	new_pred = 1 - 1 / (1+new_pred_ratio)
	#boxplot(new_pred[y_true>0.9], new_pred[y_true<0.1])
	#abline(h=0.5)
	return(new_pred)
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
