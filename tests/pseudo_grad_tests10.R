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
model.num <- 7
grphf <- models[[model.num]]

# Check the graph:
gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)
#library(Rgraphviz);dev.off()
#plot(gp)

# Adjacenty matrix:
adj <- ug(grphf, result="matrix")

s1<-1
s2<-2
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}    # Feature function

# Make up random potentials and return a CRF-object
num.samps   <- 1000
n.states    <- 2
sim.mdl    <- sim.field.random(adjacentcy.matrix=adj, num.states=n.states, num.sims=num.samps, seed=NULL)
samps       <- sim.mdl$samples
known.model <- sim.mdl$model
mrf.sample.plot(samps)


# First identify which nodes are associated with which parameters and store in the crf object:
# These are needed for the sum over k. See CRFutil for implenentation.
n2p <- nodes2params.list(known.model, storeQ = T)


psl.grads <- array(NA, c(known.model$n.par, nrow(samps)))
dim(psl.grads)
for(n in 1:nrow(samps)) {

  psl.grads[,n] <- grad.neglogpseudolik.config(
    config=samps[n,],
    phi.config=NULL,
    node.conditional.energies=NULL,
    node.complement.conditional.energies=NULL,
    par=NULL,
    crf=known.model,
    ff=f0)$gradient.log.pseudolikelihood

}
rowSums(psl.grads)
