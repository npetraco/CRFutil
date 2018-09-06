theta <- c(1.2, -3.1)
node.num <- 1
X <- c(2,2)

library(CRFutil)

# Graph formula:
grphf <- ~A:B
gp    <- ug(grphf, result = "graph")
adj   <- ug(grphf, result="matrix")
iplot(gp)
dev.off()

# Model:
n.states <- 2
fitc <- make.crf(adj, n.states)
fitc <- make.features(fitc)
fitc <- make.par(fitc, 2)
fitc$node.par[1,1,] <- 1
fitc$node.par[2,1,] <- 1
fitc$edge.par[[1]][1,1,1] <- 2
fitc$edge.par[[1]][2,2,1] <- 2

theta.reformated <- par2logpots(theta, fitc)
theta.reformated

# fitc$node.pot
# theta.reformated[[1]]
#
# fitc$node.pot
# fitc$edge.pot
# make.pots(parms=theta, crf=fitc, rescaleQ=F, replaceQ=T, printQ=F)
# fitc$node.pot
# fitc$edge.pot
# make.par.from.potentials(crf = fitc)

# Check conditional energy:
s1<-1
s2<-2
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}    # Feature function

conditional.config.energy(
  config                   = X,
  condition.element.number = node.num,
  adj.node.list            = fitc$adj.nodes,
  edge.mat                 = fitc$edges,
  one.lgp                  = theta.reformated[[1]],
  two.lgp                  = theta.reformated[[2]],
  ff                       = f0,
  printQ                   = T)


Eone(X[1], theta.reformated[[1]][[1]], f0)
Etwo(X[1], X[2], theta.reformated[[2]][[1]], f0)

#-------------------------------------------------------
#Compute phi for X_1:
phi.X <- phi.features(
  config    = X,
  edges.mat = fitc$edges,
  node.par  = fitc$node.par,
  edge.par  = fitc$edge.par,
  ff        = f0
)
phi.X

# Grab the parameters associated with conditional node (node.num)
n2p <- nodes2params.list(fitc, storeQ = T)
node.pars <- n2p[[node.num]]
node.pars

# Compute the energy: E(X_i|{\bf X}_i\slash X_i) using the dot product formulation
fitc$par <- theta
fitc$par
fitc$par[node.pars]

phi.X
phi.X[node.pars]

fitc$par[node.pars] %*% phi.X[node.pars]
