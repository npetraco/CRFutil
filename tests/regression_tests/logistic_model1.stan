data{
  int<lower=0> N;
  int<lower=-1,upper=1> X1[N];  //Delta-alpha vectors for each parameter
  int<lower=-1,upper=1> X2[N];
  int<lower=-1,upper=1> X3[N];
  int<lower=-1,upper=1> X4[N];
  int<lower=-1,upper=1> X5[N];
  int<lower=-1,upper=1> X6[N];
  int<lower=0,upper=1> y[N];   //Node responses
}
parameters{
  vector[6] beta;
}
transformed parameters{
  vector[N] Xbeta;
  for (i in 1:N)
    Xbeta[i] = beta[1]*X1[i] + beta[2]*X2[i] + beta[3]*X3[i] +
               beta[4]*X4[i] + beta[5]*X5[i] + beta[6]*X6[i];
}
model{
  beta ~ normal(0,100);
  y ~ bernoulli_logit(Xbeta);
}

