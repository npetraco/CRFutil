library(CRFutil)

load("/Users/npetraco/latex/papers/dust/steph_diss/CRFutil/tests/samps.RData")
grphf <- ~A:B + B:C + C:D + D:A
adj   <- ug(grphf, result="matrix")

n.states <- 2
#===================================
# To compare with fit:
fit <- make.crf(adj, n.states)
fit <- make.features(fit)
fit <- make.par(fit, 12)

# Fill the parameter index matrices with the standard parameterization
# Standard parameterization, one param per node:
for(i in 1:nrow(fit$node.par)){
  fit$node.par[i,1,1] <- i
}
fit$edge.par[[1]][1,1,1] <- 5
fit$edge.par[[1]][2,2,1] <- 6
fit$edge.par[[2]][1,1,1] <- 7
fit$edge.par[[2]][2,2,1] <- 8
fit$edge.par[[3]][1,1,1] <- 9
fit$edge.par[[3]][2,2,1] <- 10
fit$edge.par[[4]][1,1,1] <- 11
fit$edge.par[[4]][2,2,1] <- 12
#===================================

set.seed(1)
w.test <- runif(fit$n.par,-1,1)
w.test

# Before we use negloglik function:

fit$node.pot
fit$edge.pot

fit$par.stat <- mrf.stat(fit, samps)
fit$par.stat
fit$gradient
fit$nll
fit$par

make.pots(w.test,fit)

negloglik(w.test, fit, samps, infer.exact)

# After we use negloglik function:

# NOTE node/edge pots NOT updated. Are they WITHIN mrf.exact.nll ??????? **********
fit$node.pot
fit$edge.pot

fit$par.stat
fit$gradient
fit$nll
fit$par

infer.exact(fit)

# Update potentials with w.test and see what happens with infer.exact
upd.pots <- make.pots(parms = w.test, crf = fit, rescaleQ = F)
upd.pots[[1]]
upd.pots[[2]]

fit$node.pot <- upd.pots[[1]]
fit$edge.pot <- upd.pots[[2]]

infer.exact(fit)

# Re-run negloglik with the fit object now:
fit$par
fit$nll
negloglik(w.test, fit, samps, infer.exact)
w.test
fit$par
fit$nll
