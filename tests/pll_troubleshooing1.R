library(CRFutil)


grphf <- ~A:B
adj   <- ug(grphf, result="matrix")
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}
n.states <- 2

fitc <- make.crf(adj, n.states)
fitc <- make.features(fitc)
fitc <- make.par(fitc, 2)
fitc$node.par[1,1,] <- 1
fitc$node.par[2,1,] <- 1
fitc$node.par

fitc$edge.par[[1]][1,1,1] <- 2
fitc$edge.par[[1]][2,2,1] <- 2
fitc$edge.par

fitc$node.par
fitc$edge.par

# Choose a theta:
theta <- c(-2.8, 6.3)
theta.mats <- par2logpots(theta, fitc)
theta.mats

# Choose a config:
Xc <- c(0,1)

# Compute E(Xi|X/Xi):
conditional.config.energy(
  config                   = Xc,
  condition.element.number = 1,
  adj.node.list            = fitc$adj.nodes,
  edge.mat                 = fitc$edges,
  one.lgp                  = theta.mats[[1]],
  two.lgp                  = theta.mats[[2]],
  ff                       = f0,
  printQ                   = T)


# Checks:
make.pots(parms = theta, crf = fitc, rescaleQ = F, replaceQ = F, printQ = F)
Eone(0, theta.mats[[1]][[1]], f0)
Etwo(1, 1, theta.mats[[2]][[1]], f0)

fitc$edges
fitc$node.par
fitc$edge.par


# Compute phi(X)
Xc <- c(0,1)
phf <- phi.features(
  config    = Xc,
  edges.mat = fitc$edges,
  node.par  = fitc$node.par,
  edge.par  = fitc$edge.par,
  ff        = f0)

n2p <- nodes2params.list(fitc, storeQ = F)

cond.node.num <- 1
theta[n2p[cond.node.num][[1]]] %*% phf[n2p[cond.node.num][[1]]]
