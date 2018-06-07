# Initalize an mrf-object:
library(CRFutil)
library(Rgraphviz)


# Graph formula:
grphf <- ~A:B + B:C + C:D

# Check the graph:
gp <- ug(grphf, result = "graph")
dev.off()
plot(gp)
dev.off()
iplot(gp)

adj <- ug(grphf, result="matrix")
adj
n.states <- 2

mc <- make.crf(adj, n.states)
# These are what CRF takes as inputs/fits. The "potentials".
Psi1 <- c(0.25, 0.75)*4
Psi2 <- c(0.9,  0.1) *10
Psi3 <- c(0.25, 0.75)*4
Psi4 <- c(0.9,  0.1) *10

Psi12 <-
  6*rbind(c(2/6, 1/6),
          c(1/6, 2/6))
Psi23 <-
  6*rbind(c(2/6, 1/6),
          c(1/6, 2/6))
Psi34 <-
  6*rbind(c(2/6, 1/6),
          c(1/6, 2/6))


mc$node.pot[1,] <- Psi1
mc$node.pot[2,] <- Psi2
mc$node.pot[3,] <- Psi3
mc$node.pot[4,] <- Psi4

mc$edges # Check!
mc$edge.pot[[1]] <- Psi12
mc$edge.pot[[2]] <- Psi23
mc$edge.pot[[3]] <- Psi34

# Check again!
mc$node.pot
mc$edge.pot

samps <- sample.exact(mc, 100)
samps

# Now try to compute "sufficient statistics" and neg-log-lik
# Using standard parameterization:
mrf.fit <- make.crf(adj, n.states)
mrf.fit <- make.features(mrf.fit)
mrf.fit <- make.par(mrf.fit, 5)

mrf.fit$node.par[1,1,1]    <- 1
mrf.fit$node.par[2,1,1]    <- 2
mrf.fit$node.par[3,1,1]    <- 3
mrf.fit$node.par[4,1,1]    <- 4

for(i in 1:mrf.fit$n.edges){
  mrf.fit$edge.par[[i]][1,1,1] <- 5
  mrf.fit$edge.par[[i]][2,2,1] <- 5
}

mrf.fit$par
mrf.fit$node.pot
mrf.fit$edge.pot


# Suffient statistics:
mrf.stat(mrf.fit, samps)
UGM_MRF_computeSuffStat(mrf.fit,samps)

#
mrf.fit$node.par
mrf.fit$edge.par
UGM_MRF_makePotentials(mrf.fit$par, mrf.fit)
UGM_MRF_makePotentials(log(c(1,2,3,4,5)), mrf.fit)
mrf.fit$edge.pot
#
sfs <- UGM_MRF_computeSuffStat(mrf.fit,samps)
sfs
UGM_MRF_NLL(w = mrf.fit$par, nInstances = 100, suffStat = sfs, mrf.fit, inferFunc = infer.exact)
UGM_MRF_NLL(w = mrf.fit$par, nInstances = 100, suffStat = sfs, mrf.fit, inferFunc = infer.junction)
UGM_MRF_NLL(w = mrf.fit$par, nInstances = 100, suffStat = sfs, mrf.fit, inferFunc = infer.chain)
UGM_MRF_NLL(w = mrf.fit$par, nInstances = 100, suffStat = sfs, mrf.fit, inferFunc = infer.lbp)
#
# NLL:
mrf.fit$par.stat <- UGM_MRF_computeSuffStat(mrf.fit,samps)
mrf.nll(mrf.fit$par, mrf.fit, samps, infer.method=infer.exact)
mrf.nll(mrf.fit$par, mrf.fit, samps, infer.method=infer.junction)
#
#Gradient next

#

# Fit pots
#mrf.fit <- train.mrf(mrf.fit, nll = mrf.exact.nll, samps, infer.method = infer.exact)
#mrf.fit$nll
mrf.fit$par
mrf.fit$max.state
