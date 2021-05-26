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

phi.X <- phi.features(
  config    = X,
  edges.mat = known.model$edges,
  node.par  = known.model$node.par,
  edge.par  = known.model$edge.par,
  ff        = f0
)
phi.X

phi.Xc <- phi.features(
  config    = complement.at.idx(X,node.num),
  edges.mat = known.model$edges,
  node.par  = known.model$node.par,
  edge.par  = known.model$edge.par,
  ff        = f0
)
phi.Xc


node.pars <- n2p[[node.num]]
node.pars
phi.X[node.pars]
phi.Xc[node.pars]


conditional.config.energy(config = X,
                          condition.element.number = node.num,
                          adj.node.list = known.model$adj.nodes,
                          edge.mat = known.model$edges,
                          one.lgp = pot.info$node.energies,
                          two.lgp = pot.info$edge.energies,
                          ff = f0)

conditional.config.energy(config = complement.at.idx(X, complement.index = node.num),
                          condition.element.number = node.num,
                          adj.node.list = known.model$adj.nodes,
                          edge.mat = known.model$edges,
                          one.lgp = pot.info$node.energies,
                          two.lgp = pot.info$edge.energies,
                          ff = f0)
