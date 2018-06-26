library(CRFutil)

load("/Users/npetraco/latex/papers/dust/steph_diss/CRFutil/tests/samps.RData")
grphf <- ~A:B + B:C + C:D + D:A
adj   <- ug(grphf, result="matrix")

n.states <- 2
fit <- make.crf(adj, n.states)
fit <- make.features(fit)
fit <- make.par(fit, 12)
length(fit$par)

# Fill the parameter index matrices with the standard parameterization
# Standard parameterization, one param per node:
for(i in 1:nrow(fit$node.par)){
  fit$node.par[i,1,1] <- i
}
fit$node.par

# Use the more flexible edge parameterization:
fit$edges # Check edge order first!
fit$edge.par[[1]][1,1,1] <- 5
fit$edge.par[[1]][2,2,1] <- 6
fit$edge.par[[2]][1,1,1] <- 7
fit$edge.par[[2]][2,2,1] <- 8
fit$edge.par[[3]][1,1,1] <- 9
fit$edge.par[[3]][2,2,1] <- 10
fit$edge.par[[4]][1,1,1] <- 11
fit$edge.par[[4]][2,2,1] <- 12
fit$edge.par
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

#mrf.nll <- function(par, crf, instances, infer.method = infer.chain, ...)
set.seed(1)
w.test <- runif(com$n.par,-1,1)
w.test

com$par.stat <- mrf.stat(com, samps)
com$par.stat
com$gradient
com$nll
com$par

mrf.exact.nll(w.test, com, samps,   infer.exact)

com$node.pot
com$edge.pot
make.pots(w.test, com, rescaleQ = T)[[2]]
com$edge.pot[[2]]
make.pots(w.test, com, rescaleQ = T)[[2]][[2]]

com$edge.pot[[1]] - make.pots(w.test, com, rescaleQ = T)[[2]][[1]]
com$edge.pot[[2]] - make.pots(w.test, com, rescaleQ = T)[[2]][[2]]
com$edge.pot[[3]] - make.pots(w.test, com, rescaleQ = T)[[2]][[3]]
com$edge.pot[[4]] - make.pots(w.test, com, rescaleQ = T)[[2]][[4]]

mrf.update(com) # Rescales potentials
#mrf.exact.nll(com$par, com, samps,   infer.exact)

com$par.stat
com$gradient
com$nll
com$par


#?????
logZ1 <- infer.exact(com)$logZ
logZ1
nll.test1 <- nrow(samps)*logZ1 # cf. ll. 295 Train.cpp
nll.test1
for(i in 1:com$n.par) {        # cf. ll. 296-300 Train.cpp
  nll.test1 <- nll.test1 - w.test[i] * com$par.stat[i];
}
nll.test1
#????


fit$par.stat <- mrf.stat(fit, samps)
fit$par.stat
fit$gradient
fit$nll
fit$par

infer.exact(fit)
fit$node.pot
fp <- make.pots(w.test, fit, rescaleQ = T)
fit$node.pot <- fp[[1]]
fit$edge.pot <- fp[[2]]
fit$node.pot
fit$edge.pot
infer.exact(fit)


negloglik(w.test, fit, samps,   infer.exact)
logZ1
potsfit <- make.pots(w.test,fit)
fit$node.pot <- potsfit[[1]]
fit$edge.pot <- potsfit[[2]]
fit$edge.pot
mrf.update(fit)


UGM_MRF_NLL(w.test,50,fit$par.stat,fit,infer.exact)

fit$par.stat
fit$gradient
fit$nll
fit$par



#train.mrf <- function(crf, instances, nll = mrf.nll, infer.method = infer.chain, ..., trace = 0)
#{
#  gradient <- function(par, crf, ...) { crf$gradient }
#  crf$par.stat <- mrf.stat(crf, instances)
#  solution <- stats::optim(crf$par, nll, gradient, crf, instances, infer.method, ..., method = "L-BFGS-B", control = list(trace = trace))
#  crf$par <- solution$par
#  mrf.update(crf)
#  crf
#}
