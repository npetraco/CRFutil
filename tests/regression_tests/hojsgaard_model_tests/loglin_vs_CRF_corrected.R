library(gRbase)
library(MASS)
library(gRim)
library(CRFutil)

# Observed configurations: X
data("lizardRAW")
head(lizardRAW)

# First change the state labels to 1 or 2:
levels(lizardRAW[,1])
diam <- as.numeric(lizardRAW[,1])
head(data.frame(diam, lizardRAW[,1]))

levels(lizardRAW[,2])
hgt <- as.numeric(lizardRAW[,2])
head(data.frame(hgt, lizardRAW[,2]))

levels(lizardRAW[,3])
spc <- as.numeric(lizardRAW[,3])
head(data.frame(spc, lizardRAW[,3]))

# Desired order:
#spc  is 1
#diam is 2
#hgt  is 3

# Put the data bach together in the orignina column order and convert to a contingency table
X.RAW <- cbind(spc, diam, hgt)
head(X.RAW)
colnames(X.RAW) <- c("1","2","3")
head(X.RAW)

X <- xtabs(~., data=X.RAW)
dim(X)
ftable(X)

# First fit the Hojgaard model with loglm and gRim
ll2 <- loglm(~1:2 + 1:3, data=X);
ll2
X.loglm.coefs <- coef(ll2)
X.loglm.coefs

# Look at emp relative freqs vs the fitted relative freqs:
X.fitted <- fitted(ll2)
X.Prob.fitted <- X.fitted/sum(X.fitted)
ftable(X.Prob.fitted)

# Lizard graph
grphf <- ~1:2 + 1:3
adj <- ug(grphf, result="matrix")
adj
# Check the graph:
gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)

f0 <- function(y){ as.numeric(c((y==1),(y==2)))}
n.states <- 2

# Instantiate an empty model to fit:
knm <- make.crf(adj, n.states)
knm <- make.features(knm)
knm <- make.par(knm, 5)
knm$node.par[1,1,] <- 1
knm$node.par[2,1,] <- 2
knm$node.par[3,1,] <- 3
knm$edge.par[[1]][1,1,1] <- 4
knm$edge.par[[1]][2,2,1] <- 4
knm$edge.par[[2]][1,1,1] <- 5
knm$edge.par[[2]][2,2,1] <- 5

# Fit model to samples from the known model and obtain an estimate for w:
train.mrf(crf = knm, instances = X.RAW, nll=mrf.exact.nll, infer.method = infer.exact)
knm$par
knm$edge.pot
knm$node.pot
infer.exact(knm)


pot.info.knm <- make.gRbase.potentials(knm, node.names = gp@nodes, state.nmes = c("1","2"))
gR.dist.info.knm    <- distribution.from.potentials(pot.info.knm$node.potentials, pot.info.knm$edge.potentials)
logZ.knm            <- gR.dist.info.knm$logZ
joint.dist.info.knm <- as.data.frame(as.table(gR.dist.info.knm$state.probs))
joint.dist.info.knm
sum(joint.dist.info.knm[,4])

gR.dist.info.knm$state.probs
X.Prob.fitted
joint.dist.info.knm
joint.dist.info.knm.rearr <- joint.dist.info.knm[,c(2,3,1,4)]

X.Prob.fitted.flat <- as.data.frame(as.table(X.Prob.fitted))
joint.dist.info.knm.rearr
X.Prob.fitted.flat

# Rearrange state indices to be in the same order between the two models:
rearr.idxs <- sapply(1:8,function(xx){row.match(X.Prob.fitted.flat[xx,1:3], table = joint.dist.info.knm.rearr[,1:3])})
data.frame(X.Prob.fitted.flat[,1:3], joint.dist.info.knm.rearr[rearr.idxs,1:3])

# Compare hojsgaard loglm fit probs to CRF fit probs for lizard data:
data.frame(X.Prob.fitted.flat[,1:4], joint.dist.info.knm.rearr[rearr.idxs,4])



# Use hojsgaard loglm coefs to make potentials for CRF. Compare the state probs we get to the
# direct state probs of the hojsgaard loglm model and the original CRF model
# Compare logZ obtained to hojsgaard intercept
hoj <- make.empty.field(graph.eq=grphf, parameterization.typ="standard", plotQ=T)
dump.crf(hoj)
adj
#spc is 1
#diam is 2
#hgt is 3

hoj$node.pot[1,] <- exp(X.loglm.coefs$`1`)  # node 1 pot
hoj$node.pot[2,] <- exp(X.loglm.coefs$`2`)  # node 2 pot
hoj$node.pot[3,] <- exp(X.loglm.coefs$`3`)  # node 3 pot
hoj$node.pot

hoj$edge.pot[[1]] <- exp(t(X.loglm.coefs$`1.2`))  # 1-2
hoj$edge.pot[[2]] <- exp(t(X.loglm.coefs$`1.3`))  # 1-3
hoj$edge.pot


pot.info.hoj <- make.gRbase.potentials(hoj, node.names = gp@nodes, state.nmes = c("1","2"))
pot.info.hoj

gR.dist.info.hoj    <- distribution.from.potentials(pot.info.hoj$node.potentials, pot.info.hoj$edge.potentials)
gR.dist.info.hoj

logZ.hoj            <- gR.dist.info.hoj$logZ
logZ.hoj
X.loglm.coefs$`(Intercept)`


joint.dist.info.hoj <- as.data.frame(as.table(gR.dist.info.hoj$state.probs))
joint.dist.info.hoj

data.frame(X.Prob.fitted.flat[,1:4], joint.dist.info.hoj[rearr.idxs,4], joint.dist.info.knm.rearr[rearr.idxs,4])

# Conclusion: Hosjgaard loglm coefs ARE log potentials (energies). They are similar for those fit with CRF for this example.
