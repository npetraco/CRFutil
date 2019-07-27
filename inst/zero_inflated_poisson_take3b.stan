data {
    /* Dimensions */
    int<lower=0> N;
    int<lower=0> p;
    /* Design Matrix WITH intercept*/
    matrix[N,p] X;
    /* Outcome */
    int<lower=0> y[N];
    /* Hyperparameters*/
    real<lower=0> s_theta[p];  // was s
    real<lower=0> s_beta[p];   // was s_theta
}

parameters {
    vector[p] theta;   // was beta_lam. parameters for log lambda expansion
    vector[p] beta;    // was beta_ppi. parameters for structural-0 probability expansion
}

transformed parameters {
    vector[N] log_lambda;   // was eta
    vector[N] logit_ppi;   // was eta_theta

    log_lambda = X * theta; // log lambda was eta
    logit_ppi  = X * beta;  // logit ppi was eta_theta
}

model {
    /* Prior */
    for (j in 1:p) {
        target += normal_lpdf(theta[j] | 0, s_theta[j]);
        target += normal_lpdf(beta[j]  | 0, s_beta[j]);
    }

    /* Likelihood */
    /* https://github.com/paul-buerkner/brms/blob/master/inst/chunks/fun_zero_inflated_poisson.stan */
    for (i in 1:N) {
        if (y[i] == 0) {
            /* Zero case */
            target += log_sum_exp(/* Structural zero */
                                  bernoulli_logit_lpmf(1 | logit_ppi[i]),
                                  /* Poisson zero */
                                  bernoulli_logit_lpmf(0 | logit_ppi[i]) +
                                  poisson_log_lpmf(0 | log_lambda[i]));
        } else {
            /* Non-zero case */
            /* First term means not structural zero. */
            target += bernoulli_logit_lpmf(0 | logit_ppi[i]) +
                /* y[i] is relevant only here. */
                poisson_log_lpmf(y[i] | log_lambda[i]);
        }
    }
}
