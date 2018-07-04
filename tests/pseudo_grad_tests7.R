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


dE.mat   <- array(NA, c(known.model$n.par, known.model$n.nodes))
dE.mat.c <- array(NA, c(known.model$n.par, known.model$n.nodes))
ce.vec   <- array(NA, known.model$n.nodes)
cce.vec  <- array(NA, known.model$n.nodes)
for(i in 1:known.model$n.nodes) {

  node.num <- i
  node.pars <- n2p[[node.num]]

  dEX1.i   <- numeric(known.model$n.par)
  dEX1.i.c <- numeric(known.model$n.par)

  dEX1.i[node.pars]   <- phi.X[node.pars]
  dEX1.i.c[node.pars] <- phi.Xc[node.pars]

  dE.mat[,i]   <- dEX1.i
  dE.mat.c[,i] <- dEX1.i.c

  ce.vec[i] <- conditional.config.energy(config = X,
                            condition.element.number = node.num,
                            adj.node.list = known.model$adj.nodes,
                            edge.mat = known.model$edges,
                            one.lgp = pot.info$node.energies,
                            two.lgp = pot.info$edge.energies,
                            ff = f0)

  cce.vec[i] <- conditional.config.energy(config = complement.at.idx(X, complement.index = node.num),
                            condition.element.number = node.num,
                            adj.node.list = known.model$adj.nodes,
                            edge.mat = known.model$edges,
                            one.lgp = pot.info$node.energies,
                            two.lgp = pot.info$edge.energies,
                            ff = f0)

}
colnames(dE.mat) <- 1:known.model$n.nodes
rownames(dE.mat) <- 1:known.model$n.par
dE.mat

colnames(dE.mat.c) <- 1:known.model$n.nodes
rownames(dE.mat.c) <- 1:known.model$n.par
dE.mat.c

ce.vec
cce.vec
