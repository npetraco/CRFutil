library(CRFutil)

# Using the same node potentials (0.8,0.2), (0.5,0.5), what does the re-scaled edge potential look like under independence??
num.samps <- 100000

samps <- cbind(
  1+rbinom(n = num.samps, size = 1, prob = 0.2),
  1+rbinom(n = num.samps, size = 1, prob = 0.5)
)


mrf.sample.plot(samps)
samps

# FIT Graph formula:
grphf            <- ~A:B
adj              <- ug(grphf, result="matrix")
n.nodes        <- 2
n.states       <- 2

mrf.new <-make.crf(adj, n.states)
mrf.new <-make.features(mrf.new)
mrf.new <-make.par(mrf.new,3)

mrf.new$node.par[1,1,1] <-1
mrf.new$node.par[2,1,1] <-2
mrf.new$node.par

mrf.new$edge.par[[1]][1,1,1] <-3
mrf.new$edge.par[[1]][2,2,1] <-3
mrf.new$edge.par

mrf.new <-train.mrf(mrf.new, samps)
mrf.new$par
mrf.new$gradient

mrf.new$edge.pot[[1]]
mrf.new$node.pot/rowSums(mrf.new$node.pot)
mrf.new$edge.pot[[1]]/rowSums(mrf.new$edge.pot[[1]])
infer.exact(mrf.new)

out.mrf.new <- make.pots(parms = mrf.new$par,  crf = mrf.new,  rescaleQ = T, replaceQ = F)
mrf.new$node.pot
out.mrf.new[[1]]

mrf.new$edge.pot
out.mrf.new[[2]]

pot.info <- make.gRbase.potentials(mrf.new, node.names = c("A","B"), state.nmes = c("1","2"))
pot.info
pot.info$edge.potentials

gR.dist.info    <- distribution.from.potentials(pot.info$node.potentials, pot.info$edge.potentials)
#logZ            <- gR.dist.info$logZ
joint.dist.info <- as.data.frame(as.table(gR.dist.info$state.probs))
joint.dist.info

colnames(samps) <- c("A","B")
X <- xtabs(~., data=samps)
X
dimnames(X)

# Hojsgaard: One estimate of the config probabilities is by the relative frquencies:
X.Prob <- X/sum(X)
joint.dist.emp.info <- as.data.frame(as.table(X.Prob))

joint.dist.emp.info
joint.dist.info
library(LaplacesDemon)
KLD(joint.dist.emp.info[,3], joint.dist.info[,3])

# under independence \Psi_{ij}^{AB} are all about the same value. Scaled, they are all about 1.

# Fit with glm. What is the significance level of the interaction coef????
X.Freq <- as.data.frame(as.table(X))

glm.fit <- glm(Freq ~ -1 + A + B + A:B, family=poisson(link=log), data=X.Freq)
summary(glm.fit)

library(gRim)
dmod.fit <- dmod(~A + B + A:B, data=X)
coef(dmod.fit)

library(MASS)
ll.fit <- loglm(Freq ~ -1 + A + B + A:B, data = X.Freq)
summary(ll.fit)

glm.inp.fit <- glm(Freq ~ -1 + A + B, family=poisson(link=log), data=X.Freq)
summary(glm.inp.fit)

