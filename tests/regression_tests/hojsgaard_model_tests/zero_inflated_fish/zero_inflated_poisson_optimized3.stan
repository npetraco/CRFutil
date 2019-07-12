data {
  int<lower=0>         N_zero;               // num obs
  int<lower=0>         N_nonzero;            // num obs
  int<lower=0>         K;                    // num predictors
  matrix[N_zero, K]    x_z;                  // zeros model matrix
  matrix[N_nonzero, K] x_nz;                 // nonzeros model matrix
  int<lower=0>         y_zero[N_zero];       // zero counts
  int<lower=0>         y_nonzero[N_nonzero]; // zero counts
}
parameters {
  vector[K] beta_theta;   // Bernoulli part
  vector[K] beta_lambda;  // Poisson part
  real      alpha_theta;  // Bernoulli part
  real      alpha_lambda; // Poisson part
}
transformed parameters {
  vector[N_zero]    theta_z;       // Bernoulli prob for a zero being an inflated zero
  vector[N_nonzero] theta_nz;      // CALLED WHAT????
  vector[N_zero]    lambda_log_z;  // Poisson rate for zeros in non-inflated part
  vector[N_nonzero] lambda_log_nz; // Poisson rate for non-zeros in non-inflated part

  // insert linear predictors a + bx
  theta_z       = inv_logit(alpha_theta + x_z * beta_theta);
  theta_nz      = inv_logit(alpha_theta + x_nz * beta_theta);
  lambda_log_z  = alpha_lambda + x_z * beta_lambda;
  lambda_log_nz = alpha_lambda + x_nz * beta_lambda;
}
model {
  beta_theta   ~ normal(0,1); // Change these to cauchy
  beta_lambda  ~ normal(0,1);
  alpha_theta  ~ normal(0,1);
  alpha_lambda ~ normal(0,1);

  // Vectorized Likelihood:
  target += N_zero * log_sum_exp(bernoulli_lpmf(1 | theta_z),
                                 bernoulli_lpmf(0 | theta_z) + poisson_log_lpmf(y_zero | lambda_log_z));

  target += N_nonzero * bernoulli_lpmf(0 | theta_nz);       //THINK THE PROBLEM IS HERE. INDEX FOR NONZERO IS NOT n
  target += poisson_log_lpmf(y_nonzero | lambda_log_nz);    // NEED A SEPARATE x, theta and loglam for N_nonzero?????

}
