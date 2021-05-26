Xf[,3]

mmj <- rbind(
  c(1,  1,  1,  1),
  c(1, -1,  1, -1),
  c(1,  1, -1, -1),
  c(1, -1, -1,  1))
glm1 <- glm(Xf[,3] ~ mmj[,1] + mmj[,2] + mmj[,3] + mmj[,4] -1, family = poisson(link="log"))
coef(glm1)
megi$glm.theta.raw

library(rstanarm)
stan_glm1 <- stan_glm(Xf[,3] ~ mmj[,1] + mmj[,2] + mmj[,3] + mmj[,4] -1,
                      family = poisson,
                      data = Xf,
                      prior = normal(0,3),
                      prior_intercept = normal(0,10),
                      chains = 4,
                      cores = 4)
coef(stan_glm1)
coef(glm1)
megi$glm.theta.raw

?stan_glm
plot(stan_glm1)
pars = "beta"
plot(stan_glm1, "hist", pars = "mmj[, 1]")
plot(stan_glm1, prob=0.95, pars = "mmj[, 1]")

plot(stan_glm1, "hist", pars = "mmj[, 4]")
plot(stan_glm1, prob=0.95, pars = "mmj[, 4]")

junk <- as.matrix(stan_glm1)
head(junk)

num.samps <- 10
samps     <- sample.exact(AB, num.samps)
mrf.sample.plot(samps)
#samps
head(samps)


megb <- marginal.edge.bayes.loglin(samps)
megb$coefficients
megb$model
megb$ses
megb$y
megb$stanfit
summary(megb)
pejunk <- as.matrix(megb)
hist(pejunk[,4], xlab="omega")

megb

pepotj <- exp(pejunk[,4])/exp(-pejunk[,4])
hist(pepotj)
abline(v=1)
median(pepotj)
sd(pepotj)
