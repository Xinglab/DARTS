#include <Rcpp.h>
using namespace Rcpp;

// Re-write core functions of bayes_darts_BF.r in Rcpp for speed up
// Zijun Zhang, 1.6.2017
// last revisited: 2.2.2017
// revised 4.1.2017: add random number seed
// revised 5.18.2017: comment out some unused variables to suppress warnings in cygwin

// revised 5.13.2017: fixed the index bug in sampler
// revised 5.18.2017: comment out some unused variables to suppress warnings in cygwin


// [[Rcpp::export()]]
double dnorm_cpp(double x, double mean, double sd)
{
	return(1/sqrt(2*M_PI*pow(sd,2))*exp(-0.5*pow(x-mean,2)/pow(sd,2)));
}

// [[Rcpp::export()]]
double rdarts_post_pr_cpp(NumericVector par, NumericMatrix data, double tau, double sigma, bool is_paired=false, double inc_eff_len=2.0, double skp_eff_len=1.0)
{
	//double inc_eff_len = 2.0;
	//double skp_eff_len = 1.0;
	int n_par = par.length();
	double mu=par(0);
	double delta=par(1);
	double invalid=-pow(10,12);
	if(mu>1 || mu<0 || mu+delta>1 || mu+delta<0) {return(invalid);}
	NumericVector mu_k(n_par-2);
	// Rcout << n_par << std::endl;
	for(int i=2; i<n_par; i++)
	{
		mu_k(i-2)=par(i);
		if(par(i)>1 || par(i)<0)
		{
			return(invalid);
		}
	}
	double lik_sum=0.0;
	// prior likelihood
	double prior_norm_const = R::pnorm(1.0,0,tau,1,0) - R::pnorm(-1.0,0,tau,1,0);
	double prior = -pow(delta,2) / (2*pow(tau,2)) -log(sqrt(2*M_PI)*tau) - log(prior_norm_const);
	lik_sum = lik_sum + prior;
	//Rcout << lik_sum << std::endl;
	// group-level likelihood
	for(int i=0; i<mu_k.length(); i++)
	{
		lik_sum = lik_sum + log(dnorm_cpp(mu_k(i), mu, sigma));
	}
	//Rcout << lik_sum << std::endl;
	// individual-level likelihood
	int Inc;
	int Skp;
	int j;
	int k;
	double psi_k;
	double p_k;
	for(int i=0; i<data.nrow(); i++)
	{
		k = data(i,0);
		j = data(i,1);
		Inc=data(i,2);
		Skp=data(i,3);
		psi_k = mu_k(k-1);
		if(j==2) {psi_k=psi_k+delta;}
		if(psi_k>0.999) {psi_k=0.999;}
		if(psi_k<0.001) {psi_k=0.001;}
		p_k = inc_eff_len*psi_k/(inc_eff_len*psi_k+skp_eff_len*(1-psi_k));
		lik_sum = lik_sum + Inc*log(p_k) + Skp*log(1-p_k);	
		//Rcout << Inc*log(p_k) + Skp*log(1-p_k) << std::endl;	
	}

	//Rcout << lik_sum << std::endl;
	return(lik_sum);
}



// [[Rcpp::export()]]
NumericMatrix rdarts_posterior_MCMC_sampler(int N, double tau, double sigma, NumericMatrix data, NumericVector init_val, double proposal_width=0.01, int burnin=1000, int thinning=5, int random_state=1, bool is_paired=false, double inc_eff_len=2.0, double skp_eff_len=1.0)
{
	// initial value is crucial for shortening the burnin period.
	NumericMatrix all_sam(N, init_val.length());
	//int n_col = all_sam.ncol();
	NumericVector par_this;
	NumericVector par_next;
	double p_current;
	double p_next;
	double p_acc;
	
	all_sam(0,_)=init_val;
	for(int i=1; i<N; i++)
	{
		par_this = all_sam(i-1,_);
		par_next = par_this + rnorm(init_val.length(), 0, proposal_width);
		p_current = rdarts_post_pr_cpp(par_this, data, tau, sigma, is_paired, inc_eff_len, skp_eff_len);
		p_next = rdarts_post_pr_cpp(par_next, data, tau, sigma, is_paired, inc_eff_len, skp_eff_len);
		p_acc = exp(p_next - p_current);
		//Rcout << par_next << std::endl;
		//Rcout << p_current << std::endl;
		//Rcout << p_next << std::endl;
		//Rcout << p_acc << std::endl;
		if(runif(1)(0)<p_acc)
		{
			all_sam(i,_)=par_next;
		} else {
			all_sam(i,_)=par_this;
		}
	}
	int n_acc_sam = (N-burnin)/thinning;
	NumericMatrix acc_sam(n_acc_sam, init_val.length());
	int acc_sam_idx = 0;
	for(int i=burnin; i<N; i=i+thinning){
		acc_sam(acc_sam_idx,_)=all_sam(i,_);
		acc_sam_idx++;
	}
	
	return(acc_sam);
}


// [[Rcpp::export()]]
double darts_likelihood_cpp(double mu, double delta, List data, bool loglik, double inc_eff_len=2.0, double skp_eff_len=1.0)
{
	double invalid;
	if(loglik) {invalid=-pow(10,8);} else {invalid=0;}
	//double inc_eff_len=2.0;
	//double skp_eff_len=1.0;
	int I1=data["I1"];
	int S1=data["S1"];
	int I2=data["I2"];
	int S2=data["S2"];
	double eps=pow(10, -2);
	double psi1=mu; double psi2=mu+delta;
	//Rcout << psi1 << std::endl;
	//Rcout << psi2 << std::endl;
	if(psi2>1 || psi2<0) {return(invalid);}
	if(psi1>1 || psi1<0) {return(invalid);}
	if(psi1>1-eps) {psi1=1-eps;}
	if(psi1<eps) {psi1=eps;}
	if(psi2>1-eps) {psi2=1-eps;}
	if(psi2<eps) {psi2=eps;}
	double p1 = inc_eff_len*psi1/(inc_eff_len*psi1+skp_eff_len*(1-psi1));
	double p2 = inc_eff_len*psi2/(inc_eff_len*psi2+skp_eff_len*(1-psi2));
	//Rcout << p1 << std::endl;
	//Rcout << p2 << std::endl;
	// How strange! If i put "NumericVector L=pow...;", it does NOT work!!!
	// Need to put as "NumericVector L; L=pow...;". wtf.
	double L;
	if(loglik) {
		L = I1*log(p1) + S1*log(1-p1) + I2*log(p2) + S2*log(1-p2);
	} else {
		L=pow(p1,I1) * pow((1-p1),S1) * pow(p2,I2) * pow((1-p2),S2);
	}
	return(L);
}


/*
`set_seed` function is adapted from here: 
http://thecoatlessprofessor.com/programming/set_rs_seed_in_rcpp_sequential_case/
*/
//' Set the RNG Seed from within Rcpp
//' 
//' Within Rcpp, one can set the R session seed without triggering
//' the CRAN rng modifier check. 
//' @param seed A \code{unsigned int} that is the seed one wishes to use. 
//' @return A set RNG scope.
//' @examples
//' set.seed(10)
//' x = rnorm(5,0,1)
//' set_seed(10)
//' y = rnorm(5,0,1)
//' all.equal(x,y, check.attributes = F)

// [[Rcpp::export]]
void set_seed(unsigned int seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);  
}

// [[Rcpp::export()]]
NumericMatrix darts_posterior_MCMC_sampler(int N, double tau, List data, NumericVector init_val, double proposal_width=0.05, int burnin=1000, int thinning=1, int random_state=1, double inc_eff_len=2.0, double skp_eff_len=1.0)
{
	//set_seed(random_state);
	NumericVector mu_vec(N);
	NumericVector delta_vec(N);
	//int I1=data["I1"];
	//int S1=data["S1"];
	//int I2=data["I2"];
	//int S2=data["S2"];
	// initial value is crucial for shortening the burnin period.
	double mu_init = init_val(0);
	//Rcout << mu_init << std::endl;
	//double mu_init = 0.5;
	double delta_init = init_val(1);
	//double delta_init = 0.1;
	//Rcout << delta_init << std::endl;
	mu_vec(0) = mu_init;
	delta_vec(0) = delta_init;
	double mu_this;
	double mu_next;
	double delta_this;
	double delta_next;
	double p_current;
	double p_next;
	double p_acc;
	
	for(int i=1; i<N; i++)
	{
		mu_this = mu_vec(i-1);	
		delta_this = delta_vec(i-1);
		mu_next = mu_this + rnorm(1,0,proposal_width)(0);
		delta_next = delta_this + rnorm(1,0,proposal_width)(0);
		//p_current = dnorm_cpp(delta_this,0,tau)*darts_likelihood_cpp(mu_this,delta_this,data);
		//p_next = dnorm_cpp(delta_next,0,tau)*darts_likelihood_cpp(mu_next,delta_next,data);
		//p_acc = p_next / p_current;
		p_current = log(dnorm_cpp(delta_this,0,tau)) + darts_likelihood_cpp(mu_this,delta_this,data,true,inc_eff_len,skp_eff_len);
		p_next = log(dnorm_cpp(delta_next,0,tau)) + darts_likelihood_cpp(mu_next,delta_next,data,true,inc_eff_len,skp_eff_len);
		//Rcout << p_current << std::endl;
		//Rcout << p_next << std::endl;
		p_acc = exp(p_next - p_current);
		//Rcout << p_acc << std::endl;
		if(runif(1)(0)<p_acc)
		{
			mu_vec(i)=mu_next;
			delta_vec(i)=delta_next;
		} else {
			mu_vec(i)=mu_this;
			delta_vec(i)=delta_this;
		}
	}
	
	int n_acc_sam = (N-burnin)/thinning;
	NumericMatrix acc_sam(n_acc_sam,2);
	int acc_sam_idx = 0;
	for(int i=burnin; i<N; i=i+thinning){
		acc_sam(acc_sam_idx,0)=mu_vec(i);
		acc_sam(acc_sam_idx,1)=delta_vec(i);
		acc_sam_idx++;
	}
	
	return(acc_sam);
}
