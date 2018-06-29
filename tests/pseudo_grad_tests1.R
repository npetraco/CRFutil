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
phi.features(
  config    = samps[1,],
  edges.mat = known.model$edges,
  node.par  = known.model$node.par,
  edge.par  = known.model$edge.par,
  ff        = f0
)
#What nodes are associated with what parameter?
known.model$edges
known.model$edge.par
known.model$node.par
known.model$par

# ACCOUNT FOR NO PARAMS IN A PAR CONTAINER
j.node.par.assoc <- rep(list(NULL), known.model$n.nodes)
for(i in 1:known.model$n.nodes) {
  #print(i)
  j.par.idxs  <- as.numeric(known.model$node.par[i,,])
  j.par.idxs  <- j.par.idxs[-which(j.par.idxs == 0)]
  j.node.par.assoc[[i]] <- c(j.node.par.assoc[[i]], j.par.idxs)
  j.node.par.assoc[[i]] <- unique(j.node.par.assoc[[i]])
}
j.node.par.assoc


for(i in 1:length(known.model$edge.par)) {
  j.node.idx1 <- known.model$edges[i,1]
  j.node.idx2 <- known.model$edges[i,2]
  j.par.idxs  <- as.numeric(known.model$edge.par[[i]])
  j.par.idxs  <- j.par.idxs[-which(j.par.idxs == 0)]
  j.node.par.assoc[[j.node.idx1]] <- c(j.node.par.assoc[[j.node.idx1]], j.par.idxs)
  j.node.par.assoc[[j.node.idx2]] <- c(j.node.par.assoc[[j.node.idx2]], j.par.idxs)
  j.node.par.assoc[[j.node.idx1]] <- unique(j.node.par.assoc[[j.node.idx1]])
  j.node.par.assoc[[j.node.idx2]] <- unique(j.node.par.assoc[[j.node.idx2]])
}
j.node.par.assoc
