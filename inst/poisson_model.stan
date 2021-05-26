data{
  int p;              // Number of columns in the model matrix = number of parameters
  int N;              // Number of configurations
  int <lower=0> y[N]; // Counts of each configuration
  matrix [N,p] Mmodl; // Model matrix for the graph
}
parameters{
  vector [p] theta;   // Regression coefs.
  real   alpha;       // alpha intercept
}
model {
  // Priors
  theta ~ cauchy(0,30);
  alpha ~ cauchy(0,10);

  target += poisson_log_lpmf(y | alpha + Mmodl*theta);
}
