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

#pot.info <- make.gRbase.potentials(known.model, node.names = gp@nodes)
s1<-1
s2<-2
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}    # Feature function

# First identify which nodes are associated with which parameters and store in the crf object:
# These are needed for the sum over k. See CRFutil for implenentation.
n2p <- nodes2params.list(known.model, storeQ = T)

# Try out the new formula on the first sampled configuration, node 3:
X <- samps[1,]
X

# Make vector of phi features for the selecteted configuration
# a.
phi.vec <- phi.features(
  config    = X,
  edges.mat = known.model$edges,
  node.par  = known.model$node.par,
  edge.par  = known.model$edge.par,
  ff        = f0
)
phi.vec

# Gradient "matrix" is #parameters by #nodes. I.E., each column
# is a gradient of a condtional energy:
grad.mat   <- array(NA, c(known.model$n.par, known.model$n.nodes))
grad.c.mat <- array(NA, c(known.model$n.par, known.model$n.nodes))
# Loop over nodes:
for(i in 1:known.model$n.nodes) {

  # Definitley derivs NOT with respect to these params are 0:
  node.pars <- known.model$nodes2pars[[i]]
  # Initalize a conditional energy gradient vector for node i to 0s:
  dEX.i     <- numeric(known.model$n.par)
  dEXc.i    <- numeric(known.model$n.par)

  # Get complement phi_i
  phi.vec.c <- phi.features(
    config    = complement.at.idx(X,i),
    edges.mat = known.model$edges,
    node.par  = known.model$node.par,
    edge.par  = known.model$edge.par,
    ff        = f0
  )

  # Any phi_i = 0 in here are also 0 derivs:
  dEX.i[node.pars]  <- phi.vec[node.pars]
  dEXc.i[node.pars] <- phi.vec.c[node.pars]
  # Store gradients column-wise:
  grad.mat[,i]      <- dEX.i
  grad.c.mat[,i]    <- dEXc.i
}
colnames(grad.mat) <- 1:known.model$n.nodes
rownames(grad.mat) <- 1:known.model$n.par
colnames(grad.c.mat) <- 1:known.model$n.nodes
rownames(grad.c.mat) <- 1:known.model$n.par

grad.mat
grad.mat[7,3]
grad.c.mat
grad.c.mat[7,3]

known.model$nodes2pars[[3]]
phi.vec
X
known.model$edge.par[[2]]
known.model$edges
ii <- 4
f0(X[known.model$edges[ii,1]]) %*% known.model$edge.par[[ii]][,,1] %*% f0(X[known.model$edges[ii,2]])


# b.
node.num <- 3
ce <- conditional.config.energy(
  config                   = X,
  condition.element.number = node.num,
  adj.node.list            = known.model$adj.nodes,
  edge.mat                 = known.model$edges,
  one.lgp                  = pot.info$node.energies,
  two.lgp                  = pot.info$edge.energies,
  ff                       = f0)
ce
grad.mat[7,3]

cce <- conditional.config.energy(
  config                   = complement.at.idx(X,node.num),
  condition.element.number = node.num,
  adj.node.list            = known.model$adj.nodes,
  edge.mat                 = known.model$edges,
  one.lgp                  = pot.info$node.energies,
  two.lgp                  = pot.info$edge.energies,
  ff                       = f0)
cce
grad.c.mat[7,3]

ce * grad.mat[7,3] + cce * grad.c.mat[7,3]
