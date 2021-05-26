library(CRFutil)

# Put together known MRF model and get a sample from it:
grphf <- ~A:B + B:C + C:A
gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)

adj   <- ug(grphf, result="matrix")
adj

n.states    <- 2
tri.modl    <- sim.field.random(adjacentcy.matrix=adj, num.states=n.states, num.sims=10000, seed=1)
samps       <- tri.modl$samples
known.model <- tri.modl$model
mrf.sample.plot(samps)


fit <- make.crf(adj, n.states)
fit <- make.features(fit)
fit$node.par
fit$edge.par

fit <- make.par(fit, 6)

# One param per node:
for(i in 1:nrow(fit$node.par)){
  fit$node.par[i,1,1] <- i
}
fit$node.par

# Edge parameterization:
fit$edges # Check edge order first!
fit$edge.par[[1]][1,1,1] <- 4
fit$edge.par[[1]][2,2,1] <- 4
fit$edge.par[[2]][1,1,1] <- 5
fit$edge.par[[2]][2,2,1] <- 5
fit$edge.par[[3]][1,1,1] <- 6
fit$edge.par[[3]][2,2,1] <- 6

fit$edge.par

# Fit model to samples from the known model and obtain an estimate for w:
train.mrf(fit, samps, nll=mrf.exact.nll, infer.method = infer.exact)
fit$par.stat
theta <- fit$par
theta

# Shift potentials to w:
node.pot.orig <- fit$node.pot
edge.pot.orig <- fit$edge.pot
shift.pots(fit)
fit$node.pot
node.pot.orig

edge.pot.orig
fit$edge.pot
exp(theta)

# Now lets get Z:
infer.info <- infer.exact(fit)
logZ <- infer.info$logZ
exp(logZ)

# Define states and feature function:
s1 <- 1
s2 <- 2
f  <- function(y){ as.numeric(c((y==s1),(y==s2))) }

# Check Z
all.config <- expand.grid(c(s1,s2),c(s1,s2),c(s1,s2))
all.config
all.config.prodPots <- sapply(1:nrow(all.config), function(xx){theta %*% phi.features(all.config[xx,],fit$edges,fit$node.par,fit$edge.par,f)})
all.config.prodPots
sum(exp(all.config.prodPots)) # Same?
exp(logZ)                     # Same?

# Compute vectors of phi-features for each sampled X:
M.samp <- compute.model.matrix(samps, edges.mat = fit$edges, node.par = fit$node.par, edge.par = fit$edge.par, ff = f)
M.samp

M.all.config <- compute.model.matrix(all.config, edges.mat = fit$edges, node.par = fit$node.par, edge.par = fit$edge.par, ff = f)
M.all.config

# \text{E}_{\hat{{\boldsymbol \theta}}}[\phi_i]
# But how do we compute this without enumerating all X????:
E.thetahat.phi  <- colSums(t(sapply(1:nrow(M.all.config), function(xx){(M.all.config[xx,] * exp(theta %*% M.all.config[xx,]))/exp(logZ) })))
E.thetahat.phi

hat.E.theta.phi <- fit$par.stat/nrow(M.samp)
hat.E.theta.phi

#colSums(t(sapply(1:nrow(M.samp), function(xx){M.samp[xx,] * ((1/exp(logZ))*(exp(theta %*% M.samp[xx,])))})))

# Check NLL:
mrf.exact.nll(theta,fit,samps)
UGM_MRF_NLL(theta,10000,UGM_MRF_computeSuffStat(fit,samps),fit,infer.exact)
#fit$par.stat/nrow(M.samp) - colSums(t(sapply(1:nrow(M.samp), function(xx){(M.samp[xx,] * exp(theta %*% M.samp[xx,]))/exp(logZ) })))

# Check sufficient statistics:
fit$par.stat
colSums(M.samp)

# Check gradient:
UGM_MRF_NLL(theta,10000,UGM_MRF_computeSuffStat(fit,samps),fit,infer.exact)
fit$gradient

