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


negloglik(w.test, fit, samps, infer.exact)

# After we use negloglik function:

fit$node.pot
fit$edge.pot

fit$par.stat
fit$gradient
fit$nll
fit$par


#===================================
# To compare with fit:
com <- make.crf(adj, n.states)
com <- make.features(com)
com <- make.par(com, 12)

# Fill the parameter index matrices with the standard parameterization
# Standard parameterization, one param per node:
for(i in 1:nrow(com$node.par)){
  com$node.par[i,1,1] <- i
}
com$edge.par[[1]][1,1,1] <- 5
com$edge.par[[1]][2,2,1] <- 6
com$edge.par[[2]][1,1,1] <- 7
com$edge.par[[2]][2,2,1] <- 8
com$edge.par[[3]][1,1,1] <- 9
com$edge.par[[3]][2,2,1] <- 10
com$edge.par[[4]][1,1,1] <- 11
com$edge.par[[4]][2,2,1] <- 12
#===================================

com$node.pot
com$edge.pot

com$par.stat <- mrf.stat(com, samps)
com$par.stat
com$gradient
com$nll
com$par

mrf.exact.nll(w.test, com, samps, infer.exact)

com$node.pot
com$edge.pot

com$par.stat
fit$par.stat
com$gradient
fit$gradient
com$nll
fit$nll
com$par
fit$par

#---------------------------------------------------------
# Now try the opt (minimize the negative log likelihood):
gradient <- function(par, crf, ...) { crf$gradient }

# First with com and mrf.XXXX.nll
com$par.stat <- mrf.stat(com, samps)
mrf.exact.nll(w.test, com, samps, infer.exact)
opt.info <- stats::optim(com$par, mrf.exact.nll, gradient, com, samps, infer.exact,
                         method = "L-BFGS-B",
                         control = list(trace = 1, REPORT=1))

# Now fit fit and negloglik
fit$par.stat <- mrf.stat(fit, samps)
negloglik(w.test, fit, samps, infer.exact)
opt.info2 <- stats::optim(fit$par, negloglik, gradient, fit, samps, infer.exact,
                         method = "L-BFGS-B",
                         control = list(trace = 1, REPORT=1))


