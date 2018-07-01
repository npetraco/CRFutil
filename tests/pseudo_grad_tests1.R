library(CRFutil)

# Put together known MRF model and get a sample from it:
grphf <- ~A:B + B:C + C:A
gp <- ug(grphf, result = "graph")

adj   <- ug(grphf, result="matrix")

# Make up some potentials and get a sample of size 100:
num.samps   <- 3
n.states    <- 2
tri.modl    <- sim.field.random(adjacentcy.matrix=adj, num.states=n.states, num.sims=num.samps, seed=NULL)
samps       <- tri.modl$samples
known.model <- tri.modl$model
mrf.sample.plot(samps)
samps

# Extract parameter vector:
known.model$par      <- make.par.from.potentials(known.model)
known.model$par
# Scale the potentials to conform with the parameter vector:
rescaled.pots        <- make.pots(known.model$par, known.model, rescaleQ=FALSE, replaceQ=FALSE, printQ=F)
log(rescaled.pots[[1]])
log(rescaled.pots[[2]][[1]])
known.model$node.pot <- rescaled.pots[[1]][,,]
known.model$edge.pot <- rescaled.pots[[2]]
known.model$node.pot
known.model$edge.pot
pot.info <- make.gRbase.potentials(known.model, node.names = gp@nodes)
pot.info$node.energies
pot.info$edge.energies

s1<-1
s2<-2
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}    # Feature function


en.X1.1 <- conditional.config.energy(config        = samps[1,],
                          condition.element.number = 1,
                          adj.node.list            = known.model$adj.nodes,
                          edge.mat                 = known.model$edges,
                          one.lgp                  = pot.info$node.energies,
                          two.lgp                  = pot.info$edge.energies,
                          ff                       = f0,
                          printQ                   = F)
en.X1.1

phi.X1 <- phi.features(
  config    = samps[1,],
  edges.mat = known.model$edges,
  node.par  = known.model$node.par,
  edge.par  = known.model$edge.par,
  ff        = f0
)
phi.X1

#What nodes are associated with what parameter?
nodes2params.list(known.model, storeQ = T)
known.model$edge.par
#params2nodes.list(known.model, storeQ = T)
#known.model$par
phi.X1
known.model$node.pot
known.model$edge.pot

#known.model$nodes2pars
#known.model$pars2nodes

#
known.model$nodes2pars[[1]]
known.model$par[ known.model$nodes2pars[[1]] ]
phi.X1[ known.model$nodes2pars[[1]] ]
# Same????
known.model$par[ known.model$nodes2pars[[1]] ] %*% phi.X1[ known.model$nodes2pars[[1]] ]
en.X1.1

#
samp.num <- 3
elem.num <- 3
en.Xn.i <- conditional.config.energy(config                    = samps[samp.num,],
                                     condition.element.number = elem.num,
                                     adj.node.list            = known.model$adj.nodes,
                                     edge.mat                 = known.model$edges,
                                     one.lgp                  = pot.info$node.energies,
                                     two.lgp                  = pot.info$edge.energies,
                                     ff                       = f0,
                                     printQ                   = F)

phi.Xn <- phi.features(
  config    = samps[samp.num,],
  edges.mat = known.model$edges,
  node.par  = known.model$node.par,
  edge.par  = known.model$edge.par,
  ff        = f0
)

en.Xn.i
known.model$par[ known.model$nodes2pars[[elem.num]] ] %*% phi.Xn[ known.model$nodes2pars[[elem.num]] ]

