library(CRFutil)

load("/Users/npetraco/latex/papers/dust/steph_diss/CRFutil/tests/samps.RData")
grphf <- ~A:B + B:C + C:D + D:A
adj   <- ug(grphf, result="matrix")

n.states <- 2
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



# Try com and mrf.lbp.nll
com$par.stat <- mrf.stat(com, samps)
com$par
com$par.stat

#mrf.lbp.nll(par = com$par, crf = com, instances = samps, infer.method = infer.lbp)
#com$nll
#com$gradient

gradient <- function(par, crf, ...) { crf$gradient }
opt.info <- stats::optim(par = com$par, fn = mrf.lbp.nll, gr = gradient, crf = com, instances = samps, infer.method = infer.lbp,
                         method = "L-BFGS-B",
                         control = list(trace = 1, REPORT=1))
grad.negloglik(com, nrow(samps), com$par.stat, inference.func = infer.lbp)
com$gradient
opt.info$convergence
opt.info$message


# Optimized pots:
com$node.pot
com$edge.pot

lbp.margials.info <- infer.lbp(com)
lbp.margials.info
# > infer.lbp(com)
# $node.bel
# [,1]      [,2]
# [1,] 0.8399742 0.1600258
# [2,] 0.2200392 0.7799608
# [3,] 0.2000846 0.7999154
# [4,] 0.8199216 0.1800784
#
# $edge.bel
# $edge.bel[[1]]
# [,1]       [,2]
# [1,] 0.1199834 0.71998992
# [2,] 0.1000592 0.05996745
#
# $edge.bel[[2]]
# [,1]       [,2]
# [1,] 0.77991742 0.06005203
# [2,] 0.04000775 0.12002281
#
# $edge.bel[[3]]
# [,1]       [,2]
# [1,] 2.000729e-01 0.01997223
# [2,] 1.261672e-05 0.77994222
#
# $edge.bel[[4]]
# [,1]       [,2]
# [1,] 0.03999868 0.16007952
# [2,] 0.77992315 0.01999865
#
#
# $logZ
# [1] 9.440546

# Reset com and try com and mrf.junction.nll
com$par.stat <- mrf.stat(com, samps)
com$par
com$par.stat

#mrf.junction.nll(par = com$par, crf = com, instances = samps, infer.method = infer.junction)
#com$nll
#com$gradient

gradient <- function(par, crf, ...) { crf$gradient }
opt.info <- stats::optim(par = com$par, fn = mrf.junction.nll, gr = gradient, crf = com, instances = samps, infer.method = infer.junction,
                         method = "L-BFGS-B",
                         control = list(trace = 1, REPORT=1))
com$gradient
opt.info$convergence
opt.info$message

# Optimized pots:
com$node.pot
com$edge.pot

junction.margials.info <-infer.junction(com)
# $node.bel
# [,1]      [,2]
# [1,] 0.8413270 0.1586730
# [2,] 0.2362173 0.7637827
# [3,] 0.2194966 0.7805034
# [4,] 0.8081201 0.1918799
#
# $edge.bel
# $edge.bel[[1]]
# [,1]       [,2]
# [1,] 0.1244947 0.71683233
# [2,] 0.1117226 0.04695039
#
# $edge.bel[[2]]
# [,1]       [,2]
# [1,] 0.76804416 0.07328287
# [2,] 0.04007594 0.11859704
#
# $edge.bel[[3]]
# [,1]       [,2]
# [1,] 0.213813589 0.02240369
# [2,] 0.005683038 0.75809968
#
# $edge.bel[[4]]
# [,1]       [,2]
# [1,] 0.04357209 0.17592454
# [2,] 0.76454801 0.01595537
#
#
# $logZ
# [1] 6.559518


# Reset com and try com and mrf.exact.nll
com$par.stat <- mrf.stat(com, samps)
com$par
com$par.stat

mrf.exact.nll(par = com$par, crf = com, instances = samps, infer.method = infer.exact)
com$nll
com$gradient

gradient <- function(par, crf, ...) { crf$gradient }
opt.info <- stats::optim(par = com$par, fn = mrf.exact.nll, gr = gradient, crf = com, instances = samps, infer.method = infer.exact,
                         method = "L-BFGS-B",
                         control = list(trace = 1, REPORT=1))
opt.info$convergence
opt.info$message

# Optimized pots:
com$node.pot
com$edge.pot

exact.margials.info <-infer.exact(com)
exact.margials.info
# > infer.exact(com)
# $node.bel
# [,1]      [,2]
# [1,] 0.8400072 0.1599928
# [2,] 0.2199997 0.7800003
# [3,] 0.1999913 0.8000087
# [4,] 0.8200066 0.1799934
#
# $edge.bel
# $edge.bel[[1]]
# [,1]       [,2]
# [1,] 0.12000183 0.72000535
# [2,] 0.09999787 0.05999496
#
# $edge.bel[[2]]
# [,1]       [,2]
# [1,] 0.78001253 0.05999465
# [2,] 0.03999404 0.11999879
#
# $edge.bel[[3]]
# [,1]      [,2]
# [1,] 1.999913e-01 0.0200084
# [2,] 1.267991e-08 0.7800003
#
# $edge.bel[[4]]
# [,1]       [,2]
# [1,] 0.03999878 0.15999253
# [2,] 0.78000778 0.02000091
#
#
# $logZ
# [1] 31.62207
