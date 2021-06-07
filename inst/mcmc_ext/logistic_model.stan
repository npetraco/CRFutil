data{
  int<lower=0> N;
  int<lower=0> K;
  matrix<lower=-1,upper=1> [N, K] Delta_alpha;
  int<lower=0,upper=1> y[N];   //Node responses
}
parameters{
  vector[K] theta;
}
model{
  theta ~ normal(0,100);
  //y ~ bernoulli_logit(Delta_alpha*theta);
  target += bernoulli_logit_lpmf(y | Delta_alpha*theta);
}
