library(rstan)

# Simulate some fake data
a.tru <- 3.2
b1.tru <- 10.76
sig.tru <- 0.12

alps <- rnorm(100, mean = a.tru,  sd = 10)
bets <- rnorm(100, mean = b1.tru, sd = 1.5)
sigs <- abs(rnorm(100, mean = sig.tru, sd = 0.01))

x   <- sample(c(-1,1), 100, replace = T)
mus <- alps + bets*x
y   <- sapply(1:100, function(xx){rnorm(1, mean = mus[xx], sd = sigs[xx])})

plot(x,y)
# Standard MLE linear regression:
mle.fit <- lm(y ~ x)
summary(mle.fit)
abline(mle.fit)


x.new   <- sample(c(-1,1), 100, replace = T)
dat <- list(N    = 100,
            x    = x,
            y    = y,
            Nnew = 100,
            xnew = x.new
)


# Compile and run the model:
stan.c <-"
data {
  int<lower=0> N;
  vector[N]    x;
  vector[N]    y;
  int<lower=0> Nnew;
  vector[Nnew] xnew;
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}
model {
  alpha ~ cauchy(0,1);
  beta  ~ cauchy(0,5);
  sigma ~ normal(0,1);
  y ~ normal(alpha + beta * x, sigma);
}
generated quantities {
  vector[Nnew] mu_pred;
  vector[Nnew] y_tilde;

  mu_pred = alpha + beta * xnew;

  for (i in 1:Nnew){
    y_tilde[i] = normal_rng(mu_pred[i], sigma);
  }
}"

#rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores())

# Compile and run the model
model.c <- stanc(model_code = stan.c, model_name = 'model1')
sm <- stan_model(stanc_ret = model.c, verbose = T)
bayes.fit <- sampling(sm, data = dat, iter=5000, thin = 1, chains = 4)
#bayes.fit <- sampling(sm, data = dat, iter=5000, thin = 1, chains = 4, control = list(adapt_delta = 0.9))
bayes.fit

# Examine the sampling output in more detail:
alpha <- extract(bayes.fit,"alpha")[[1]]
plot(alpha, typ="l", main="Trace")
hist(alpha, bre=80, probability = T)
acf(alpha)
mean(alpha)
median(alpha)
a.tru

beta <- extract(bayes.fit,"beta")[[1]]
plot(beta, typ="l", main="Trace")
hist(beta, bre=80, probability = T)
acf(beta)
mean(beta)
median(beta)
b1.tru

sigma <- extract(bayes.fit,"sigma")[[1]]
plot(sigma, typ="l", main="Trace")
hist(sigma, bre=80, probability = T)
acf(sigma)
mean(sigma)
median(sigma)
sig.tru

