data {
  int<lower=0> N;    // num obs
  int<lower=0> K;    // num predictors
  matrix[N, K] x;    // model matrix
  int<lower=0> y[N]; // counts
}
parameters {
  vector[K] beta_theta;   // Regression coefs for Bernoulli part
  vector[K] beta_lambda;  // Regression coefs Poisson part
  real      alpha_theta;  // Intercept for Bernoulli part
  real      alpha_lambda; // Intercept for Poisson part
}
transformed parameters {
  vector[N] theta;       // Bernoulli prob for a zero being an inflated zero
  vector[N] lambda_log;  // Poisson rate for non-inflated part

  theta      =  inv_logit(alpha_theta + x * beta_theta); // insert linear predictors a + bx
  lambda_log =  alpha_lambda + x * beta_lambda;
}
model {
  // beta_theta   ~ normal(0,100); // Change these to cauchy
  // beta_lambda  ~ normal(0,100);
  // alpha_theta  ~ normal(0,100);
  // alpha_lambda ~ normal(0,100);
  beta_theta   ~ cauchy(0,5); // Change these to cauchy
  beta_lambda  ~ cauchy(0,5);
  alpha_theta  ~ cauchy(0,5);
  alpha_lambda ~ cauchy(0,5);


  // Likelihood:
  // Pr(y|lambda,theta) = theta + (1-theta)*exp(-lambda)  if y == 0
  //                      (1-theta)*poisson(lambda)       if y > 0
  for(n in 1:N){
    if(y[n] == 0){ // for zeros: zero inflated or poisson part
      target += log_sum_exp(bernoulli_lpmf(1 | theta[n]),
                            bernoulli_lpmf(0 | theta[n]) + poisson_log_lpmf(y[n] | lambda_log[n]));
    } else {      // for non-zeros: poisson part
      target += bernoulli_lpmf(0 | theta[n]) + poisson_log_lpmf(y[n] | lambda_log[n]);
    }

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
	    log_lik[n] =   bernoulli_lpmf(0 | theta[n]) + poisson_log_lpmf(y[n] | lambda_log[n]);
    }
  }

}
