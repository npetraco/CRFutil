library(CRFutil)

# Put together known MRF model and get a sample from it:
grphf <- ~A:B + B:C + C:A
gp <- ug(grphf, result = "graph")

adj   <- ug(grphf, result="matrix")

# Make up some potentials and get a sample of size 100:
num.samps   <- 3
n.states    <- 2
tri.modl    <- sim.field.random(adjacentcy.matrix=adj, num.states=n.states, num.sims=num.samps, seed=6)
samps       <- tri.modl$samples
known.model <- tri.modl$model
mrf.sample.plot(samps)
samps

s1<-1
s2<-2
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}    # Feature function
pot.info <- make.gRbase.potentials(known.model, node.names = gp@nodes)


conditional.config.energy(config                   = samps[1,],
                          condition.element.number = 3,
                          adj.node.list            = known.model$adj.nodes,
                          edge.mat                 = known.model$edges,
                          one.lgp                  = pot.info$node.energies,
                          two.lgp                  = pot.info$edge.potentials,
                          ff                       = f0,
                          printQ                   = F)

known.model$adj.nodes
phi1 <- phi.features(
  config    = samps[1,],
  edges.mat = known.model$edges,
  node.par  = known.model$node.par,
  edge.par  = known.model$edge.par,
  ff        = f0
)
known.model$par
#What nodes are associated with what parameter?
nodes2params.list(known.model, storeQ = T)
params2nodes.list(known.model, storeQ = T)

