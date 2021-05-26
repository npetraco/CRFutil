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

#pot.info <- make.gRbase.potentials(known.model, node.names = gp@nodes)
s1<-1
s2<-2
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}    # Feature function

# First identify which nodes are associated with which parameters and store in the crf object:
# These are needed for the sum over k. See CRFutil for implenentation.
n2p <- nodes2params.list(known.model, storeQ = T)

# Try out the new formula on the first sampled configuration, node 3:
X <- samps[1,]
X

phi.X <- phi.features(config = t(as.matrix(X)),
                      edges.mat = known.model$edges, node.par = known.model$node.par,
                      edge.par = known.model$edge.par, ff = f0)
phi.X

conditional.config.energy2(
  config                   = t(as.matrix(X)),
  phi.config               = phi.X,
  condition.element.number = 3,
  crf                      = known.model,
  ff                       = f0)

known.model$par

conditional.config.energy(
  config                   = X,
  condition.element.number = 3,
  adj.node.list            = known.model$adj.nodes,
  edge.mat                 = known.model$edges,
  one.lgp                  = pot.info$node.energies,
  two.lgp                  = pot.info$edge.energies,
  ff                       = f0)

# ???????????????
neglogpseudolik.config(config = X, crf = known.model, ff = f0)

psl.info$condtional.energies
psl.info$complement.condtional.energies
psl.info$pseudo.likelihoods
psl.info$conditional.Prs
psl.info$complement.conditional.Prs

-log(psl.info$pseudo.likelihoods)
