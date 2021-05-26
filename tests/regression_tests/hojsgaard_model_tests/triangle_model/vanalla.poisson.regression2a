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
  real <lower=0> lambda[N];

  for (i in 1:N) {
    // Linear predictor
    lp[i] = dot_product(beta,Xmodl[i]);

    // Mean
    lambda[i] = exp(lp[i]);
  }
}
model {
  // Prior on regression coefs
  beta ~ cauchy(0,30);

  // Likelihood part of Bayesian inference
  y ~ poisson(lambda);
}
// generated quantities{
//   vector[N] logZ;
//   for(i in 1:N) {
//   	logZ[i] = log(N) - log(lambda[i]) + lp[i];
//   }
// }
