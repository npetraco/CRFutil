library(CRFutil)

# Put together known MRF model and get a sample from it:
grphf <- ~A:B + B:C + C:A
gp <- ug(grphf, result = "graph")

adj   <- ug(grphf, result="matrix")

# Make up some potentials and get a sample of size 100:
num.samps   <- 100
n.states    <- 2
tri.modl    <- sim.field.random(adjacentcy.matrix=adj, num.states=n.states, num.sims=num.samps, seed=NULL)
samps       <- tri.modl$samples
known.model <- tri.modl$model
mrf.sample.plot(samps)

# Extract parameter vector:
known.model$par      <- make.par.from.potentials(known.model)

# Scale the potentials to conform with the parameter vector:
rescaled.pots        <- make.pots(known.model$par, known.model, rescaleQ=FALSE, replaceQ=T, printQ=F)
pot.info <- make.gRbase.potentials(known.model, node.names = gp@nodes)


s1<-1
s2<-2
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}    # Feature function


#What nodes are associated with what parameter?
nodes2params.list(known.model, storeQ = T)

#
en.result <- NULL
for(i in 1: num.samps) {
  samp.num <- i

  phi.Xn <- phi.features(
    config    = samps[samp.num,],
    edges.mat = known.model$edges,
    node.par  = known.model$node.par,
    edge.par  = known.model$edge.par,
    ff        = f0
  )

  for(j in 1:ncol(samps)) {
    elem.num <- j

    en.Xn.i <- conditional.config.energy(config                    = samps[samp.num,],
                                         condition.element.number = elem.num,
                                         adj.node.list            = known.model$adj.nodes,
                                         edge.mat                 = known.model$edges,
                                         one.lgp                  = pot.info$node.energies,
                                         two.lgp                  = pot.info$edge.energies,
                                         ff                       = f0,
                                         printQ                   = F)
    en.phi.i.j <- known.model$par[ known.model$nodes2pars[[elem.num]] ] %*% phi.Xn[ known.model$nodes2pars[[elem.num]] ]
    #print(paste("samp#=",i,"ele#=",j,"Formula.1=",en.Xn.i,"Formula.2=",en.phi.i.j, "Delta=", en.Xn.i - en.phi.i.j))
    en.result <- rbind(en.result, c(i,j,en.Xn.i,en.phi.i.j, en.Xn.i - en.phi.i.j))

  }
}
en.result
en.result[,5]

samp.num <- 2
elem.num <- 3
conditional.config.energy2(config = samps[samp.num,], condition.element.number = elem.num, crf = known.model, ff = f0)
conditional.config.energy(config= samps[samp.num,],condition.element.number = elem.num,adj.node.list= known.model$adj.nodes,edge.mat= known.model$edges,one.lgp= pot.info$node.energies,two.lgp= pot.info$edge.energies,ff= f0,printQ= F)

