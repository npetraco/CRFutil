library(CRFutil)

# Graph formulas:
models <- list(
  ~A:B,
  ~A:B + B:C + C:D,
  ~A:B:C,
  ~A:B:C:D,
  ~A:B + B:C + C:D + D:A,
  ~A:B:C + C:D,
  ~A:B + A:C + A:D + A:E + B:C + B:D + B:E + C:D + D:E,
  ~A:B:C + C:X + D:E:G + G:X + Q:R:S + S:X,
  ~A:B + B:C + C:D + E:F + F:G + G:H + A:E + B:F + C:G + D:H + I:J + J:K + K:L + I:E + J:F + K:G + L:H + M:N + N:O + O:P + M:I + N:J + O:K + L:P)


# Choose a model
model.num <- 4
grphf <- models[[model.num]]

# Check the graph:
gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)
#library(Rgraphviz);dev.off()
#plot(gp)

# Adjacenty matrix:
adj <- ug(grphf, result="matrix")


# Make up random potentials and return a CRF-object
num.samps   <- 100
n.states    <- 2
sim.mdl    <- sim.field.random(adjacentcy.matrix=adj, num.states=n.states, num.sims=num.samps, seed=NULL)
samps       <- sim.mdl$samples
known.model <- sim.mdl$model
mrf.sample.plot(samps)

pot.info <- make.gRbase.potentials(known.model, node.names = gp@nodes)
s1<-1
s2<-2
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}    # Feature function

# First identify which nodes are associated with which parameters and store in the crf object:
# These are needed for the sum over k. See CRFutil for implenentation.
n2p <- nodes2params.list(known.model, storeQ = T)

# Try out the new formula on the first sampled configuration, node 3:
#X <- samps[1,]
X <- samps[sample(1:num.samps, size = 1),]
X

# Make vector of phi features for the selecteted configuration
# a.
phi.vec <- phi.features(
  config    = X,
  edges.mat = known.model$edges,
  node.par  = known.model$node.par,
  edge.par  = known.model$edge.par,
  ff        = f0
)
phi.vec

# Gradient "matrix" is #parameters by #nodes. I.E., each column
# is a gradient of a condtional energy:
grad.Ec.mat   <- array(NA, c(known.model$n.par, known.model$n.nodes))
grad.Ec.c.mat <- array(NA, c(known.model$n.par, known.model$n.nodes))
# Condtional and complement condtional energy:
ce            <- array(NA, c(known.model$n.nodes))
cce           <- array(NA, c(known.model$n.nodes))
# Condtional partition functions:
Zc            <- array(NA, c(known.model$n.nodes))
# Gradients of condtional partition functions:
dZc.mat       <- array(NA, c(known.model$n.par, known.model$n.nodes))

# Compute all parts of pseudo-likelihood gradient:
# Loop over nodes:
for(i in 1:known.model$n.nodes) {

  # **Definitley derivs NOT with respect to these params are 0:
  node.pars <- known.model$nodes2pars[[i]]
  # Initalize a conditional energy gradient vectors for node i to 0s:
  dEX.i     <- numeric(known.model$n.par)
  dEXc.i    <- numeric(known.model$n.par)

  # Get complement phi_i's needed for Zcs and dZcs:
  phi.vec.c <- phi.features(
    config    = complement.at.idx(X,i),
    edges.mat = known.model$edges,
    node.par  = known.model$node.par,
    edge.par  = known.model$edge.par,
    ff        = f0
  )

  # **Any phi_i = 0 in here are also 0 derivs:
  dEX.i[node.pars]  <- phi.vec[node.pars]
  dEXc.i[node.pars] <- phi.vec.c[node.pars]
  # Store condional energy gradients column-wise:
  grad.Ec.mat[,i]      <- dEX.i
  grad.Ec.c.mat[,i]    <- dEXc.i

  ce[i] <- conditional.config.energy(
    config                   = X,
    condition.element.number = i,
    adj.node.list            = known.model$adj.nodes,
    edge.mat                 = known.model$edges,
    one.lgp                  = pot.info$node.energies,
    two.lgp                  = pot.info$edge.energies,
    ff                       = f0)

  cce[i] <- conditional.config.energy(
    config                   = complement.at.idx(X,i),
    condition.element.number = i,
    adj.node.list            = known.model$adj.nodes,
    edge.mat                 = known.model$edges,
    one.lgp                  = pot.info$node.energies,
    two.lgp                  = pot.info$edge.energies,
    ff                       = f0)

  # Store conditional Z gradients column-wise as well:
  # NOTE: Assumes only two states per node
  dZc.mat[, i] <- exp(ce[i]) * grad.Ec.mat[,i] + exp(cce[i]) * grad.Ec.c.mat[,i]

  # Condtional Zs:
  # NOTE: Assumes only two states per node
  Zc[i]        <- exp(ce[i]) + exp(cce[i])

}


# GENERALIZE THIS TO WORK WITH conditional.config.energy2
psl.info <- pseudolikelihoods.from.energies(
  t(as.matrix(X)),
  known.model$adj.nodes,
  known.model$edges,
  pot.info$node.energies,
  pot.info$edge.energies,
  conditional.config.energy,
  f0)

# What happens with these if we have all possible configs??
# Average features with respect to states of X_i
E.Xphi.mat <- array(NA, c(known.model$n.par, known.model$n.nodes))
for(i in 1:known.model$n.nodes) {
  E.Xphi.mat[,i] <- dZc.mat[,i]/Zc[i]
}
E.Xphi.mat

# psl grad for a config????
grad.psl.X <- rowSums(grad.Ec.mat - E.Xphi.mat)
grad.psl.X


#??????????????
junk <- grad.neglogpseudolik.config(config=X,
                            phi.config=NULL,
                            node.conditional.energies=NULL,
                            node.complement.conditional.energies=NULL,
                            par=NULL,
                            crf=known.model,
                            ff=f0)

junk$conditional.energies
ce
junk$complement.conditional.energies
cce

junk$gradients.conditional.energies - grad.Ec.mat
junk$gradients.complement.conditional.energies - grad.Ec.c.mat

junk$gradients.conditional.partition.functions
dZc.mat
junk$gradients.conditional.partition.functions - dZc.mat

junk$gradients.log.conditional.partition.functions
E.Xphi.mat
junk$gradients.log.conditional.partition.functions - E.Xphi.mat

junk$gradient.log.pseudolikelihood
grad.psl.X
junk$gradient.log.pseudolikelihood - grad.psl.X

