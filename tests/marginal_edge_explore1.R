#library(CRF)
library(CRFutil)
#library(gRbase)

n.nodes        <- 2
n.states       <- 2
prior.prob     <- c(0.8,0.2)
trans.prob     <- matrix(0,nrow=2,ncol=2)
trans.prob[1,] <- c(0.95,0.05)
trans.prob[2,] <- c(0.05,0.95)

# FIT Graph formula:
grphf            <- ~A:B
adj              <- ug(grphf, result="matrix")
mc               <- make.crf(adj, n.states)
mc$node.pot[1,]  <- prior.prob
#mc$node.pot[2,]  <- prior.prob
mc$edge.pot[[1]] <- trans.prob

num.samps <- 100000
samps     <- sample.exact(mc, num.samps)
mrf.sample.plot(samps)
samps

mrf.new <-make.crf(adj, n.states)
mrf.new <-make.features(mrf.new)
mrf.new <-make.par(mrf.new,3)

mrf.new$node.par[1,1,1] <-1
mrf.new$node.par[2,1,1] <-2
mrf.new$node.par

mrf.new$edge.par[[1]][1,1,1] <-3
mrf.new$edge.par[[1]][2,2,1] <-3
#mrf.new$edge.par[[1]][1,2,1] <-3
#mrf.new$edge.par[[1]][2,1,1] <-4

mrf.new$edge.par

mrf.new <-train.mrf(mrf.new, samps)

mrf.new$par
mrf.new$node.pot/rowSums(mrf.new$node.pot)
mrf.new$edge.pot[[1]]/rowSums(mrf.new$edge.pot[[1]])

infer.exact(mrf.new)

# Now get dmod estimates

# Now get glm estimates

# True config probs
mc.true <-make.crf(adj, n.states)
mc.true <-make.features(mc.true)
mc.true <-make.par(mc.true,3)

mc.true$node.par[1,1,1] <-1
mc.true$node.par[2,1,1] <-2
mc.true$node.par

mc.true$edge.par[[1]][1,1,1] <- 3
mc.true$edge.par[[1]][2,2,1] <- 3
mc.true$edge.par

mc.true$node.pot[1,]  <- prior.prob
#mc.true$node.pot[2,]  <- prior.prob
mc.true$edge.pot[[1]] <- trans.prob

theta.true <- make.par.from.potentials(mc.true)
theta.true

out.mc <- make.pots(parms = theta.true,  crf = mc.true,  rescaleQ = T, replaceQ = T)
out.mc
mc.true$par <- theta.true
mc.true$node.pot
mc.true$edge.pot

pot.true.info <- make.gRbase.potentials(mc.true, node.names = c("A","B"), state.nmes = c("1","2"))
pot.true.info

gR.dist.true.info    <- distribution.from.potentials(pot.true.info$node.potentials, pot.true.info$edge.potentials)
#logZ            <- gR.dist.info$logZ
joint.dist.true.info <- as.data.frame(as.table(gR.dist.true.info$state.probs))
joint.dist.true.info

# Empirical configuration probs from sample:
# Fold the observed configurations into a contingency table
colnames(samps) <- c("A","B")
X <- xtabs(~., data=samps)
X
dimnames(X)

# Hojsgaard: One estimate of the config probabilities is by the relative frquencies:
X.Prob <- X/sum(X)
X.Prob
joint.dist.emp.info <- as.data.frame(as.table(X.Prob))
joint.dist.emp.info



# Check FIT config probs
out.mrf.new <- make.pots(parms = mrf.new$par,  crf = mrf.new,  rescaleQ = T, replaceQ = T)
mrf.new$node.pot
mrf.new$edge.pot
mrf.new$par

pot.info <- make.gRbase.potentials(mrf.new, node.names = c("A","B"), state.nmes = c("1","2"))
pot.info
pot.info$edge.potentials

gR.dist.info    <- distribution.from.potentials(pot.info$node.potentials, pot.info$edge.potentials)
#logZ            <- gR.dist.info$logZ
joint.dist.info <- as.data.frame(as.table(gR.dist.info$state.probs))
joint.dist.info

# from bels get cond with gR

# use kl dist between emp dist and fit dist to assess fit between models??
