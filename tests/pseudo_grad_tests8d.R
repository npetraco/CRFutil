library(CRFutil)

# Graph formula for Slayer field:
grphf <- ~A:B+A:C+A:D+A:E+B:C+B:D+B:E+C:D+D:E

# Check the graph:
gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)

# Adjacenty matrix:
adj <- ug(grphf, result="matrix")

#pot.info <- make.gRbase.potentials(known.model, node.names = gp@nodes)
s1<-1
s2<-2
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}    # Feature function

# Make up random potentials and return a CRF-object
num.samps   <- 100
n.states    <- 2
slay    <- sim.field.random(adjacentcy.matrix=adj, num.states=n.states, num.sims=num.samps, seed=NULL)
samps       <- slay$samples
known.model <- slay$model
mrf.sample.plot(samps)

# First identify which nodes are associated with which parameters and store in the crf object:
# These are needed for the sum over k. See CRFutil for implenentation.
n2p <- nodes2params.list(known.model, storeQ = T)

# Compute pseudo-likelihoods and intermediates using featue function based formulation:
pot.info <- make.gRbase.potentials(known.model, node.names = gp@nodes)
psl.info <- pseudolikelihoods.from.energies(samps,known.model$adj.nodes,known.model$edges,pot.info$node.energies,
  pot.info$edge.energies,conditional.config.energy,f0)
psl.info


# neg log pseudolikelihoods for each sampled config using featue function based formulation:
-log(psl.info$pseudo.likelihoods)
# from feature formulation
nlpslks <-  sapply(1:nrow(samps), function(xx){neglogpseudolik.config(config = samps[xx,], crf = known.model, ff = f0)})
nlpslks
-log(psl.info$pseudo.likelihoods) - nlpslks
dev.off()
plot(-log(psl.info$pseudo.likelihoods) - nlpslks, typ="h")
