library(CRFutil)

# Graph formula for Slayer field:
grphf <- ~A:B
gp    <- ug(grphf, result = "graph")
adj   <- ug(grphf, result="matrix")

# Make up random potentials and return a CRF-object
# num.samps   <- 100
# n.states    <- 2
# slay        <- sim.field.random(adjacentcy.matrix=adj, num.states=n.states,
#                                 num.sims=num.samps, seed=1)
# samps       <- slay$samples
# known.model <- slay$model
# mrf.sample.plot(samps)
n.states <- 2
fitc <- make.crf(adj, n.states)
fitc <- make.features(fitc)
fitc <- make.par(fitc, 2)
fitc$node.par[1,1,] <- 1
fitc$node.par[2,1,] <- 1
fitc$edge.par[[1]][1,1,1] <- 2
fitc$edge.par[[1]][2,2,1] <- 2


# Choose a theta:
theta <- c(-2.8, 6.3)
theta.reformated <- par2logpots(theta, fitc)
fitc$par <- theta

#pot.info <- make.gRbase.potentials(known.model, node.names = gp@nodes)
s1<-1
s2<-2
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}    # Feature function

# First identify which nodes are associated with which parameters and store in the crf object:
# These are needed for the sum over k. See CRFutil for implementation.
n2p <- nodes2params.list(fitc, storeQ = T)
fitc$nodes2pars

# Try out the new formula on the first sampled configuration, node 3:
X <- c(2,1)
node.num <- 1

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
node.pars <- n2p[[node.num]]
node.pars

# Compute the energy: E(X_i|{\bf X}_i\slash X_i) using the dot product formulation
fitc$par[node.pars] %*% phi.X[node.pars]

# Compute the energy: E(X_i|{\bf X}_i\slash X_i) using the origninal recipie (mmmmmm... KFC)
conditional.config.energy(
  config                   = X,
  condition.element.number = node.num,
  adj.node.list            = fitc$adj.nodes,
  edge.mat                 = fitc$edges,
  one.lgp                  = theta.reformated[[1]],
  two.lgp                  = theta.reformated[[2]],
  ff                       = f0,
  printQ                   = T)


conditional.config.energy(config=X,
                           condition.element.number = node.num,
                           adj.node.list= fitc$adj.nodes,
                           edge.mat= fitc$edges,
                           one.lgp= theta.reformated[[1]],
                           two.lgp= theta.reformated[[2]],
                           ff= f0,printQ= T)



