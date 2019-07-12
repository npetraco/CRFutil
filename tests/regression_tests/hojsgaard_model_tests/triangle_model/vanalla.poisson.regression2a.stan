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
transformed parameters {

  real lp[N];
  real <lower=0> lambda[N];

  for (i in 1:N) {
    // Linear predictor
    lp[i] = alpha + dot_product(theta,Mmodl[i]);

    // Mean
    lambda[i] = exp(lp[i]);
  }
}
model {
  // Prior on regression coefs
  theta ~ cauchy(0,30);
  alpha ~ cauchy(0,10);

  // Likelihood part of Bayesian inference
  y ~ poisson(lambda);
}
