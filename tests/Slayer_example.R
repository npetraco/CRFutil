library(CRFutil)

# Graph formula for Slayer field:
grphf <- ~A:B+A:C+A:D+A:E+B:C+B:D+B:E+C:D+D:E

# Check the graph:
gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)

# Adjacenty matrix:
adj <- ug(grphf, result="matrix")

# Define states and feature function:
s1 <- 1
s2 <- 2
f  <- function(y){ as.numeric(c((y==s1),(y==s2))) }

# Make up random potentials and return a CRF-object
slay <- sim.field.random(adjacentcy.matrix=adj, num.states=2, num.sims=1, seed=1)$model

# Get energies from potentials anf decorate both with gRBase annotations:
slay.params   <- make.gRbase.potentials(crf=slay, node.names=gp@nodes, state.nmes=c(s1,s2))
slay.node.en  <- slay.params$node.energies
slay.edge.en  <- slay.params$edge.energies
slay.node.pot <- slay.params$node.potentials
slay.edge.pot <- slay.params$edge.potentials

# Enumerate all the state configurations
config.mat <- expand.grid(c(s1,s2),c(s1,s2),c(s1,s2),c(s1,s2),c(s1,s2))
colnames(config.mat) <- c("A","B","C","D","E")

# Compute Probs of all configurations and logZ from energies:
dist.en.info <- distribution.from.energies(
  state.space = config.mat,
  edges.mat = slay$edges,
  node.energies = slay.node.en,
  edge.energies = slay.edge.en,
  energy.func = config.energy,
  ff = f)
Pr.en <- dist.en.info$state.probs

# As a check compute Probs of all configurations and logZ from potentials as well:
dist.pot.info <- distribution.from.potentials(
  gRbase.node.potentials = slay.node.pot,
  gRbase.edge.potentials = slay.edge.pot)
Pr.pot <- dist.pot.info$state.probs

# Check quick to make sure all config energies are the same between the two combination methods:
library(prodlim)
names(dimnames(Pr.pot))
order(names(dimnames(Pr.pot)))

Pr.pot2 <- as.data.frame(as.table(Pr.pot))
Pr.pot2
gR.idxs <- Pr.pot2[,order(names(dimnames(Pr.pot)))]

# Rearrange gR base state order to be in the same order as config.mat:
rearrange.idxs <- sapply(1:nrow(config.mat), function(xx){row.match(config.mat[xx,], table = gR.idxs)})
rearrange.idxs

# Columns the same?
cbind(Pr.pot2[rearrange.idxs,6],Pr.en, Pr.pot2[rearrange.idxs,6]-Pr.en)

# Put in a nice table:
prodPots        <- Pr.en*exp(dist.en.info$logZ) # work backwards
config.energies <- log(prodPots)                # work backwards
cbind(config.mat, config.energies, prodPots, Pr.en, Pr.pot2[rearrange.idxs,6])
