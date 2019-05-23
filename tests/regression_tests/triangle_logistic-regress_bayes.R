library(CRFutil)
library(rstan)
library(shinystan)
library(coda)

# Fully connected
grphf <- ~1:2+1:3+2:3
adj <- ug(grphf, result="matrix")
# Check the graph:
gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)

f0 <- function(y){ as.numeric(c((y==1),(y==2)))}
n.states <- 2

# Instantiate the "true model":
true.model <- make.crf(adj, n.states)
true.model <- make.features(true.model)
true.model <- make.par(true.model, 6)
true.model$node.par[1,1,] <- 1
true.model$node.par[2,1,] <- 2
true.model$node.par[3,1,] <- 3
true.model$edge.par[[1]][1,1,1] <- 4
true.model$edge.par[[1]][2,2,1] <- 4
true.model$edge.par[[2]][1,1,1] <- 5
true.model$edge.par[[2]][2,2,1] <- 5
true.model$edge.par[[3]][1,1,1] <- 6
true.model$edge.par[[3]][2,2,1] <- 6

set.seed(6)
true.model$par <- runif(6,-1.5,1.1)
true.model$par # "true" theta
out.pot <- make.pots(parms = true.model$par,  crf = true.model,  rescaleQ = T, replaceQ = T)


# So now sample from the model as if we obtained an experimental sample:
num.samps <- 25
set.seed(1)
samps <- sample.exact(true.model, num.samps)
mrf.sample.plot(samps)

# Instantiate an empty model for the Bayesian fit:
bayes.lrm <- make.crf(adj, n.states)
bayes.lrm <- make.features(bayes.lrm)
bayes.lrm <- make.par(bayes.lrm, 6)
bayes.lrm$node.par[1,1,] <- 1
bayes.lrm$node.par[2,1,] <- 2
bayes.lrm$node.par[3,1,] <- 3

bayes.lrm$edge.par[[1]][1,1,1] <- 4
bayes.lrm$edge.par[[1]][2,2,1] <- 4
bayes.lrm$edge.par[[2]][1,1,1] <- 5
bayes.lrm$edge.par[[2]][2,2,1] <- 5
bayes.lrm$edge.par[[3]][1,1,1] <- 6
bayes.lrm$edge.par[[3]][2,2,1] <- 6
bayes.lrm$edges

#----------
mle.lrm <- make.crf(adj, n.states)
mle.lrm <- make.features(mle.lrm)
mle.lrm <- make.par(mle.lrm, 6)
mle.lrm$node.par[1,1,] <- 1
mle.lrm$node.par[2,1,] <- 2
mle.lrm$node.par[3,1,] <- 3

mle.lrm$edge.par[[1]][1,1,1] <- 4
mle.lrm$edge.par[[1]][2,2,1] <- 4
mle.lrm$edge.par[[2]][1,1,1] <- 5
mle.lrm$edge.par[[2]][2,2,1] <- 5
mle.lrm$edge.par[[3]][1,1,1] <- 6
mle.lrm$edge.par[[3]][2,2,1] <- 6
#----------

Delta.alpha.info <- delta.alpha(crf = bayes.lrm, samples = samps, printQ = F)
Delta.alpha <- Delta.alpha.info$Delta.alpha

y <-c(samps[,1], samps[,2], samps[,3])
y
y[which(y==2)] <- 0
y

dat <- list(
  N=length(y),
  X1=Delta.alpha[,1],
  X2=Delta.alpha[,2],
  X3=Delta.alpha[,3],
  X4=Delta.alpha[,4],
  X5=Delta.alpha[,5],
  X6=Delta.alpha[,6],
  y = y
)

model.c <- stanc(file =
            "/home/npetraco/codes/R/CRFutil/tests/regression_tests/logistic_model1.stan",
            model_name = 'model')
sm <- stan_model(stanc_ret = model.c, verbose = T)

# Run the model
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
bfit <- sampling(sm, data = dat,
                #control = list(adapt_delta = 0.85),
                iter=20000,
                thin = 10,
                chains = 8)
print(bfit)


beta <- extract(bfit,"beta")[[1]]
dim(beta)
# plot(beta[,1], typ="l", main="Trace")
# hist(beta[,1], bre=80, probability = T)
# acf(beta[,1])
# mean(beta[,1])
# median(beta[,1])

launch_shinystan(bfit)

Ma <- glm(y ~ Delta.alpha[,1] + Delta.alpha[,2] + Delta.alpha[,3] +
              Delta.alpha[,4] + Delta.alpha[,5] + Delta.alpha[,6] -1,
              family=binomial(link="logit"))
summary(Ma)


# c(mean(beta[,1]), mean(beta[,2]), mean(beta[,3]), mean(beta[,4]), mean(beta[,5]), mean(beta[,6]))
# c(median(beta[,1]), median(beta[,2]), median(beta[,3]), median(beta[,4]), median(beta[,5]), median(beta[,6]))
# true.model$par

bayes.lrm$par <- c(median(beta[,1]), median(beta[,2]), median(beta[,3]), median(beta[,4]), median(beta[,5]), median(beta[,6]))
mle.lrm$par <- as.numeric(coef(Ma))
mle.lrm$par

out.pot2 <- make.pots(parms = bayes.lrm$par,  crf = bayes.lrm,  rescaleQ = T, replaceQ = T)
bayes.lrm$node.pot
bayes.lrm$edge.pot

out.pot2a <- make.pots(parms = mle.lrm$par,  crf = mle.lrm,  rescaleQ = T, replaceQ = T)
mle.lrm$node.pot
mle.lrm$edge.pot


# Node and edge beliefs:
bayes.lrm.bel  <- infer.exact(bayes.lrm)
mle.lrm.bel <- infer.exact(mle.lrm)
true.model.bel  <- infer.exact(true.model)

bayes.lrm.bel$node.bel
mle.lrm.bel$node.bel
true.model.bel$node.bel

bayes.lrm.bel$edge.bel[[1]]
mle.lrm.bel$edge.bel[[1]]
true.model.bel$edge.bel[[1]]


bayes.lrm.bel$edge.bel[[2]]
mle.lrm.bel$edge.bel[[2]]
true.model.bel$edge.bel[[2]]


bayes.lrm.bel$edge.bel[[3]]
mle.lrm.bel$edge.bel[[3]]
true.model.bel$edge.bel[[3]]

# Configuration probabilities:
pot.info <- make.gRbase.potentials(bayes.lrm, node.names = gp@nodes, state.nmes = c("1","2"))
pot.info

gR.dist.info    <- distribution.from.potentials(pot.info$node.potentials, pot.info$edge.potentials)
logZ            <- gR.dist.info$logZ
joint.dist.info <- as.data.frame(as.table(gR.dist.info$state.probs))
joint.dist.info

pot.info.true.model <- make.gRbase.potentials(true.model, node.names = gp@nodes, state.nmes = c("1","2"))
pot.info.true.model


pot.info2 <- make.gRbase.potentials(mle.lrm, node.names = gp@nodes, state.nmes = c("1","2"))
pot.info2

gR.dist.info2    <- distribution.from.potentials(pot.info2$node.potentials, pot.info2$edge.potentials)
logZ2            <- gR.dist.info2$logZ
joint.dist.info2 <- as.data.frame(as.table(gR.dist.info2$state.probs))
joint.dist.info2

pot.info.true.model <- make.gRbase.potentials(true.model, node.names = gp@nodes, state.nmes = c("1","2"))


gR.dist.info.true.model    <- distribution.from.potentials(pot.info.true.model$node.potentials, pot.info.true.model$edge.potentials)
logZ.true.model            <- gR.dist.info.true.model$logZ
joint.dist.info.true.model <- as.data.frame(as.table(gR.dist.info.true.model$state.probs))
joint.dist.info.true.model

bayes.lrm.cp  <- round(joint.dist.info[,4]*100, 1)     # Bayes logistic
mle.lrm.cp <- round(joint.dist.info2[,4]*100, 1)    # MLE logistic
true.model.cp  <- round(joint.dist.info.true.model[,4]*100, 1) # True
cbind(joint.dist.info[,c(2,3,1)], bayes.lrm.cp, mle.lrm.cp, true.model.cp)

sum(bayes.lrm.cp)
sum(true.model.cp)

process.logistic.fit.stan(bfit)
junk <- HPDinterval(as.mcmc(beta), prob = 0.95)
class(junk)
junk

#mrf.nll(par = mle.lrm$par, crf = mle.lrm, instances = t(as.matrix(c(1,1,1))), infer.method = infer.exact)
mle.lrm.bel$logZ
logZ2


# Some checks:

# Enumerate all the state configurations
s1 <- 1
s2 <- 2
all.configs <- expand.grid(c(s1,s2),c(s1,s2),c(s1,s2))
colnames(all.configs) <- gp@nodes
all.configs

fit.params.info <- make.gRbase.potentials(mle.lrm, node.names=gp@nodes, state.nmes=c(s1,s2))
fit.params.info


config.energies <- sapply(1:nrow(all.configs),
                          function(xx){
                            config.energy(config = all.configs[xx,], edges.mat = mle.lrm$edges,
                                          one.lgp = fit.params.info$node.energies,
                                          two.lgp = fit.params.info$edge.energies,
                                          ff = f0)
                          })
config.energies
logZ2

# Compare config energies with energies computed as theta \dot phi:
# First compute the “features” (phi) for all each possible configuration
M.all  <- compute.model.matrix(all.configs, mle.lrm$edges, mle.lrm$node.par, mle.lrm$edge.par, f0)
M.all

# Now get the energy of each config with the dot-porduct formula:
w <- mle.lrm$par
alt.config.energies <- sapply(1:nrow(M.all), function(xx){w%*%M.all[xx,]})
config.energies                       # Check: All config energies with feature functions
alt.config.energies                   # Check: All config energies with features


#   ?????? did we re-normalize at some point?????

mle.lrma <- copy.crf(mle.lrm, plotQ = F)

out.pot3 <- make.pots(parms = mle.lrma$par,  crf = mle.lrma,  rescaleQ = F, replaceQ = T)
lz.unsc <- infer.exact(crf = mle.lrma)$logZ
alt.config.energies - lz.unsc
config.energies - logZ2

# yup

config.energies - alt.config.energies # Check: Differences


# So be careful when computing config probs that you are using the correct normalization constant!
# IE did you re-scale the potentials?????

conditional.config.energy()
