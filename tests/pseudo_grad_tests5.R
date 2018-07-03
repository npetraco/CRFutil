library(CRFutil)

# Graph formula for Slayer field:
grphf <- ~A:B+A:C+A:D+A:E+B:C+B:D+B:E+C:D+D:E

# Check the graph:
gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)

# Adjacenty matrix:
adj <- ug(grphf, result="matrix")


# Make up random potentials and return a CRF-object
num.samps   <- 100
n.states    <- 2
slay    <- sim.field.random(adjacentcy.matrix=adj, num.states=n.states, num.sims=num.samps, seed=1)
samps       <- slay$samples
known.model <- slay$model
mrf.sample.plot(samps)

pot.info <- make.gRbase.potentials(known.model, node.names = gp@nodes)
s1<-1
s2<-2
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}    # Feature function

# First identify which nodes are associated with which parameters and store in the crf object:
# These are needed for the sum over k. See CRFutil for implenentation.
n2p <- nodes2params.list(known.model, storeQ = T)

# Try out the new formula on the first sampled configuration, node 3:
X <- samps[1,]
node.num <- 3

#Compute phi for X_1:
phi.X <- phi.features(
  config    = X,
  edges.mat = known.model$edges,
  node.par  = known.model$node.par,
  edge.par  = known.model$edge.par,
  ff        = f0
)
phi.X


# Grab the parameters associated with node 3
node.pars <- n2p[[node.num]]
node.pars        # ** Definitley derivs not with respect to these params are 0
phi.X[node.pars] # ** Any phi_i = 0 in here are also 0 derivs

# Compute the energy:
# E(X_1|{\bf X}\slash X_1) = \theta_1 \phi_1({\bf X}) + \sum_{ k \in \{3,7,10,13\} } \theta_k \phi_{k_{[1\sim j\in\{2,3,4,5\}]}}({\bf X})
known.model$par[node.pars] %*% phi.X[node.pars]

# Gradient at this energy: \frac{\partial}{\partial \theta_k} E(X_i|{\bf X}\slash X_i) = \phi_{k_{[ \sim i]}}({\bf X})
# Gradient at this energy ????:
dEX1.3 <- numeric(known.model$n.par)
dEX1.3[node.pars] <- phi.X[node.pars]
dEX1.3
node.pars
phi.X[node.pars]

dE.mat <- NULL
for(i in 1:known.model$n.nodes) {

  node.num <- i
  node.pars <- n2p[[node.num]]
  dEX1.i <- numeric(known.model$n.par)
  dEX1.i[node.pars] <- phi.X[node.pars]
  print(dEX1.i)
  dE.mat <- cbind(dE.mat, dEX1.i)
}
colnames(dE.mat) <- 1:known.model$n.nodes
rownames(dE.mat) <- 1:known.model$n.par
dE.mat
# CHECK??
