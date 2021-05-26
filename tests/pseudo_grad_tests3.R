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

pot.info <- make.gRbase.potentials(known.model, node.names = gp@nodes)
s1<-1
s2<-2
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}    # Feature function

en.result <- array(0, c(num.samps*ncol(samps),5))
count <- 1
for(i in 1:num.samps) {
  for(j in 1:ncol(samps)) {
    samp.num <- i
    elem.num <- j
    ce1 <- conditional.config.energy2(config = samps[samp.num,], condition.element.number = elem.num, crf = known.model, ff = f0)
    ce2 <- conditional.config.energy(config= samps[samp.num,],condition.element.number = elem.num,adj.node.list= known.model$adj.nodes,edge.mat= known.model$edges,one.lgp= pot.info$node.energies,two.lgp= pot.info$edge.energies,ff= f0,printQ= F)
    en.result[count,] <- c(i,j,ce1,ce2,ce1-ce2)
    count <- count + 1
  }
}
en.result
en.result[,5]
