functions {
  int num_zeros(int[] y) {
    int ssum = 0;
    for (n in 1:size(y))
      ssum += (y[n] == 0);
    return ssum;
  }
}
data {
  int<lower=0> N;    // num obs
  int<lower=0> K;    // num predictors
  matrix[N, K] x;    // model matrix
  int<lower=0> y[N]; // counts
}
transformed data {
  int<lower = 0> N_zero = num_zeros(y);
  int<lower = 1> y_nonzero[N - N_zero];
  int N_nonzero = 0;
  for (n in 1:N) {
    if (y[n] == 0) continue;
    N_nonzero += 1;
    y_nonzero[N_nonzero] = y[n];
  }
}
parameters {
  vector[K] beta_theta;   // Bernoulli part
  vector[K] beta_lambda;  // Poisson part
  real      alpha_theta;  // Bernoulli part
  real      alpha_lambda; // Poisson part
}
transformed parameters {
  vector[N] theta;       // Bernoulli prob for a zero being an inflated zero
  vector[N] lambda_log;  // Poisson rate for non-inflated part

  theta      =  inv_logit(alpha_theta + x * beta_theta); // insert linear predictors a + bx
  lambda_log =  alpha_lambda + x * beta_lambda;
}
model {
  beta_theta   ~ normal(0,1); // Change these to cauchy
  beta_lambda  ~ normal(0,1);
  alpha_theta  ~ normal(0,1);
  alpha_lambda ~ normal(0,1);

  // Likelihood:
  target += N_zero * log_sum_exp(bernoulli_lpmf(1 | theta),
                                 bernoulli_lpmf(0 | theta) + poisson_log_lpmf(y | lambda_log));

  target += N_nonzero * bernoulli_lpmf(0 | theta);       //THINK THE PROBLEM IS HERE. INDEX FOR NONZERO IS NOT n
  target += poisson_log_lpmf(y_nonzero | lambda_log);

}
