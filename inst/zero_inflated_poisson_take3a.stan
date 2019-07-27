data {
    /* Dimensions */
    int<lower=0> N;
    int<lower=0> M;
    /* Design Matrix */
    matrix[N,M] X;
    /* Outcome */
    int<lower=0> y[N];
    /* Hyperparameters*/
    real<lower=0> s[M];
    real<lower=0> s_theta[M];
}

parameters {
    vector[M] beta;
    vector[M] beta_theta;
}

transformed parameters {
    vector[N] eta;
    vector[N] eta_theta;

    eta = X * beta;
    eta_theta = X * beta_theta;
}

model {
    /* Prior */
    for (j in 1:M) {
        target += normal_lpdf(beta[j] | 0, s[j]);
        target += normal_lpdf(beta_theta[j] | 0, s_theta[j]);
    }

    /* Likelihood */
    /* https://github.com/paul-buerkner/brms/blob/master/inst/chunks/fun_zero_inflated_poisson.stan */
    for (i in 1:N) {
        if (y[i] == 0) {
            /* Zero case */
            target += log_sum_exp(/* Structural zero */
                                  bernoulli_logit_lpmf(1 | eta_theta[i]),
                                  /* Poisson zero */
                                  bernoulli_logit_lpmf(0 | eta_theta[i]) +
                                  poisson_log_lpmf(0 | eta[i]));
        } else {
            /* Non-zero case */
            /* First term means not structural zero. */
            target += bernoulli_logit_lpmf(0 | eta_theta[i]) +
                /* y[i] is relevant only here. */
                poisson_log_lpmf(y[i] | eta[i]);
        }
    }
}
