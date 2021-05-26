library(CRFutil)

# Put together known MRF model and get a sample from it:
grphf <- ~A:B + B:C + C:A
gp <- ug(grphf, result = "graph")

adj   <- ug(grphf, result="matrix")

# Make up some potentials and get a sample of size 100:
num.samps   <- 100
n.states    <- 2
tri.modl    <- sim.field.random(adjacentcy.matrix=adj, num.states=n.states, num.sims=num.samps, seed=1)
samps       <- tri.modl$samples
known.model <- tri.modl$model
mrf.sample.plot(samps)

# Using the sample, fit a model in the standard parameterization to obtain a theta:
fit <- make.crf(adj, n.states)
fit <- make.features(fit)
fit$node.par # Node parameter index matrix (empty)
fit$edge.par # Edge parameter index matrices (empty)

fit <- make.par(fit, 6)

# Fill the parameter index matrices with the standard parameterization
# Standard parameterization, one param per node:
for(i in 1:nrow(fit$node.par)){
  fit$node.par[i,1,1] <- i
}
fit$node.par # Check

# Standard parameterization, edge parameterization:
fit$edges # Check edge order first!
fit$edge.par[[1]][1,1,1] <- 4
fit$edge.par[[1]][2,2,1] <- 4
fit$edge.par[[2]][1,1,1] <- 5
fit$edge.par[[2]][2,2,1] <- 5
fit$edge.par[[3]][1,1,1] <- 6
fit$edge.par[[3]][2,2,1] <- 6

fit$edge.par # Check

#------------------------------------------------
train.mrf(fit, samps, nll=mrf.exact.nll, infer.method = infer.exact)
# train.mrf re-scales the potentials as a last setp. Put them back the way they were:
shift.pots(fit) # Necessary??

#fit$node.pot # probably shifted by mrf.update
#fit$node.par # there is no indication of this in node.par however

# Compute marginals:
infr.info <- infer.exact(fit)
node.bels <- infr.info$node.bel
edge.bels <- infr.info$edge.bel

node.phi.means <- array(0,fit$n.nodes)
for(i in 1:nrow(fit$node.par)){
  for(j in 1:ncol(fit$node.par)) {
    print(paste(infr.info$node.bel[i,j] , fit$node.par[i,j,1], (1 - fit$node.par[i,j,1]<=0)))
    node.phi.means[i]  <- node.phi.means[i]  + infr.info$node.bel[i,j] * (1 - fit$node.par[i,j,1]<=0)
  }
}
node.phi.means
infr.info$node.bel

edge.phi.means <- array(0,fit$n.edges)
for(k in 1:fit$n.edges){
  for(i in 1:nrow(fit$edge.par[[k]])){
    for(j in 1:ncol(fit$edge.par[[k]])) {

      print(paste(fit$edge.par[[k]][i,j,1], edge.bels[[k]][i,j]))
      edge.phi.means[k] <- edge.phi.means[k] + edge.bels[[k]][i,j] * (1 - fit$edge.par[[k]][i,j,1]<=0)
    }
  }
}
edge.phi.means

c(node.phi.means, edge.phi.means)
feature.means(fit, inference.func = infer.exact) # Check

# from orig pots
#0.2400000 0.4100000 0.6600001 0.4700000 0.3000000 0.2500000
# from update.pot scaled pots:
#0.2399992 0.4100003 0.6599998 0.4700005 0.3000011 0.2499983
