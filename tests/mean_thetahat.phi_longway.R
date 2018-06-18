library(CRFutil)

# Put together known MRF model and get a sample from it:
grphf <- ~A:B + B:C + C:A
gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)

adj   <- ug(grphf, result="matrix")
adj

# Make up some potentials and get a sample of size 100:
num.samps   <- 100
n.states    <- 2
tri.modl    <- sim.field.random(adjacentcy.matrix=adj, num.states=n.states, num.sims=num.samps, seed=1)
samps       <- tri.modl$samples
known.model <- tri.modl$model
mrf.sample.plot(samps)

# Using the sample, fit a model in the standard parameterization to obtain a theta:
fit <- make.crf(adj, n.states)
fit <- make.features(fit)
fit$node.par # Node parameter index matrix (empty)
fit$edge.par # Edge parameter index matrices (empty)

fit <- make.par(fit, 6)

# Fill the parameter index matrices with the standard parameterization
# Standard parameterization, one param per node:
for(i in 1:nrow(fit$node.par)){
  fit$node.par[i,1,1] <- i
}
fit$node.par

# Standard parameterization, edge parameterization:
fit$edges # Check edge order first!
fit$edge.par[[1]][1,1,1] <- 4
fit$edge.par[[1]][2,2,1] <- 4
fit$edge.par[[2]][1,1,1] <- 5
fit$edge.par[[2]][2,2,1] <- 5
fit$edge.par[[3]][1,1,1] <- 6
fit$edge.par[[3]][2,2,1] <- 6

fit$edge.par

# Fit model to samples from the known model and obtain an estimate for theta:
train.mrf(fit, samps, nll=mrf.exact.nll, infer.method = infer.exact)
fit$par.stat     # Sample sufficient statistics: total number of appearances of each parameter in the sample
theta <- fit$par # \hat{theta}
theta

# mrf.update within train.mrf switches the order of the parameters.
# Shift potentials back to be in the same order as in node.par and edge.par.
# If they are not in the same order prodPots will differ depending on the formula used to compute them.
# Only necessary for testing purposes:
shift.pots(fit)

# Compute marginals and Z (inference):
infr.info <- infer.exact(fit)

# Define states and feature function:
s1 <- 1
s2 <- 2
f  <- function(y){ as.numeric(c((y==s1),(y==s2))) }

# Enumerate all possible configurations X and compute their phi vectors:
# Note: Generally you don't want to do this for larger networks...
X.all <- expand.grid(c(s1,s2),c(s1,s2),c(s1,s2))
M.all <- compute.model.matrix(X.all, fit$edges, fit$node.par, fit$edge.par, f)

# Compute exp(E(X)) for all X as exp(\hat{theta} \dot phi(X)):
prodPots <- sapply(1:nrow(M.all), function(xx){exp(theta %*% M.all[xx,])})
prodPots

#Compute config probs:
Prs <- (1/exp(infr.info$logZ)) * prodPots
Prs
sum(Prs)

# Finally compute E_{\hat{\theta}_i}[phi_i(X)] = \sum_X \phi_i(X) Pr(X)
E.hattheta.phi <- colSums(sapply(1:ncol(M.all), function(xx){M.all[,xx] * Prs})) # Estimate at the optimal theta
hatE.theta.phi <- fit$par.stat/num.samps                                         # Empirical estimate: \hat{E}_{\theta_i}[\phi_i]
grad           <- hatE.theta.phi - E.hattheta.phi                                # gradient at the optimum theta
cbind(E.hattheta.phi,hatE.theta.phi,grad)
