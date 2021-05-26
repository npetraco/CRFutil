data{
  int p;              // Number of columns in the model matrix = number of parameters
  int N;              // Number of configurations
  int <lower=0> y[N]; // Counts of each configuration
  real Mmodl[N,p];    // Model matrix for the graph
}
parameters{
  real theta[p];       // Regression coefs.
  real alpha;          // alpha
}
model {
  // Prior on regression coefs
  theta ~ cauchy(0,30);
  alpha ~ cauchy(0,10);

  // Likelihood
  for (i in 1:N) {
    //y[i] ~ poisson_log(alpha + dot_product(theta,Mmodl[i]));
    target += poisson_log_lpmf(y[i] | alpha + dot_product(theta,Mmodl[i]));
  }
  //y ~ poisson_log(alpha + Mmodl*theta); // Why doesn't this work

}
