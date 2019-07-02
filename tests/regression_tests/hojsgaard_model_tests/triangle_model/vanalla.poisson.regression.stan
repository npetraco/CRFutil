data{
  int p;              // Number of columns in the model matrix
  int N;              // Number of observed configurations
  int <lower=0> y[N]; // Counts of each observed config
  real Xmodl[N,p];    // Model matrix for the graph
}
parameters{
  real beta[p];       // Regression coefs.
}
transformed parameters {

  real lp[N];
  real <lower=0> mu[N];

  for (i in 1:N) {
    // Linear predictor
    lp[i] = dot_product(beta,Xmodl[i]);

    // Mean
    mu[i] = exp(lp[i]);
  }
}
model {
  // Prior part of Bayesian inference
  // Flat prior for mu (no need to specify if non-informative)
  beta ~ normal(0,30);

  // Likelihood part of Bayesian inference
  y ~ poisson(mu);
}
