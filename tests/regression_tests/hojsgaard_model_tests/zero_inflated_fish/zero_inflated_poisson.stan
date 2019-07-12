data {
  int<lower=0> N;
  int<lower=0> K; 
  matrix[N, K] x; 
  int<lower=0> y[N];
}
parameters {
  vector[K] beta_theta;                 //
  vector[K] beta_lambda;                 //
  real      alpha_theta;
  real      alpha_lambda;
}
transformed parameters {
  vector[N] theta;
  vector[N] lambda_log;
  
  theta      =  inv_logit(alpha_theta + x * beta_theta); 
  lambda_log =  alpha_lambda + x * beta_lambda; 
}
model {
  beta_theta   ~ normal(0,1); 
  beta_lambda  ~ normal(0,1); 
  alpha_theta  ~ normal(0,1);
  alpha_lambda ~ normal(0,1);
  
  for(n in 1:N){
    if(y[n] == 0)
      target += log_sum_exp(bernoulli_lpmf(1 | theta[n]),
                            bernoulli_lpmf(0 | theta[n]) + poisson_log_lpmf(y[n] | lambda_log[n]));
    else
      target += bernoulli_lpmf(0 | theta[n])
                  + poisson_log_lpmf(y[n] | lambda_log[n]);
  }
}
generated quantities{
  real log_lik[N];
  int<lower=0> y_sim[N];
  int zero;
  for(n in 1:N){
    zero = bernoulli_rng(theta[n]);
	  y_sim[n] = (1-zero)*poisson_log_rng(lambda_log[n]);
	     
	  if(y[n] == 0){
	    log_lik[n] = log_sum_exp(bernoulli_lpmf(1 | theta[n]),
                               bernoulli_lpmf(0 | theta[n]) + poisson_log_lpmf(y[n] | lambda_log[n])) ;
	  } else {
	    log_lik[n] =   bernoulli_lpmf(0 | theta[n])
                     + poisson_log_lpmf(y[n] | lambda_log[n]);
    }
  }
  
}
