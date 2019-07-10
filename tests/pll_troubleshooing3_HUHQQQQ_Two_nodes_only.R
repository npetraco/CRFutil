library(CRFutil)

# Graph formula for Slayer field:
grphf <- ~A:B
gp    <- ug(grphf, result = "graph")
adj   <- ug(grphf, result="matrix")
iplot(gp)
dev.off()

# Make up random potentials and return a CRF-object
num.samps   <- 100
n.states    <- 2
slay        <- sim.field.random.simplest(adjacentcy.matrix=adj, num.states=n.states,
                                num.sims=num.samps, seed=1)
samps       <- slay$samples
known.model <- slay$model
mrf.sample.plot(samps)

# Needed for the energy functions:
pot.info <- make.gRbase.potentials(known.model, node.names = gp@nodes)
s1<-1
s2<-2
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}    # Feature function

# First identify which nodes are associated with which parameters and store in the crf object:
# These are needed for the sum over k. See CRFutil for implementation.
n2p <- nodes2params.list(known.model, storeQ = T)

# Try out the new formula on the first sampled configuration, node 3:
X <- samps[1,]
node.num <- 2

#Compute phi for X_1:
phi.X <- phi.features(
  config    = X,
  edges.mat = known.model$edges,
  node.par  = known.model$node.par,
  edge.par  = known.model$edge.par,
  ff        = f0
)

# Grab the parameters associated with node 3
node.pars <- n2p[[node.num]]
node.pars

# Compute the energy: E(X_1|{\bf X}_1\slash X_3)
known.model$par[node.pars] %*% phi.X[node.pars]

# Compute the energy with the original rcipie (mmmmm KFC!):
conditional.config.energy(config=X,
                          condition.element.number = node.num,
                          adj.node.list= known.model$adj.nodes,
                          edge.mat= known.model$edges,
                          one.lgp= pot.info$node.energies,
                          two.lgp= pot.info$edge.energies,
                          ff= f0,printQ= F)

# For storage:
en.result <- array(0, c(num.samps*ncol(samps),5))

# Loop around the sample numbers and then the elements of each sample
count <- 1
for(i in 1:num.samps) {
  for(j in 1:ncol(samps)) {
    samp.num <- i
    elem.num <- j
    # Feature formulation:
    ce1 <- conditional.config.energy2(config = samps[samp.num,],
                                      condition.element.number = elem.num, crf =
                                        known.model, ff = f0)
    # Feature function formulation:
    ce2 <- conditional.config.energy(config=samps[samp.num,],
                                     condition.element.number = elem.num,
                                     adj.node.list= known.model$adj.nodes,
                                     edge.mat= known.model$edges,
                                     one.lgp= pot.info$node.energies,
                                     two.lgp= pot.info$edge.energies,
                                     ff= f0,printQ= F)

    en.result[count,] <- c(i,j,ce1,ce2,ce1-ce2)
    count <- count + 1
  }
}
en.result
en.result[,5]
plot(1:length(en.result[,5]),en.result[,5], ylab="Difference", xlab="Config index")
