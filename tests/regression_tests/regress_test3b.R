library(CRFutil)

# Fully connected
grphf <- ~1:2+1:3+2:3
adj <- ug(grphf, result="matrix")
# Check the graph:
gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)

f0 <- function(y){ as.numeric(c((y==1),(y==2)))}
n.states <- 2

# Instantiate an empty model to fit:
knm <- make.crf(adj, n.states)
knm <- make.features(knm)
knm <- make.par(knm, 6)
knm$node.par[1,1,] <- 1
knm$node.par[2,1,] <- 2
knm$node.par[3,1,] <- 3
knm$edge.par[[1]][1,1,1] <- 4
knm$edge.par[[1]][2,2,1] <- 4
knm$edge.par[[2]][1,1,1] <- 5
knm$edge.par[[2]][2,2,1] <- 5
knm$edge.par[[3]][1,1,1] <- 6
knm$edge.par[[3]][2,2,1] <- 6

#set.seed(6)
knm$par <- runif(6,-1.5,1.1)
knm$par # "true" theta
out.pot <- make.pots(parms = knm$par,  crf = knm,  rescaleQ = T, replaceQ = T)

# So now sample from the model as if we obtained an experimental sample:
num.samps <- 25
#set.seed(1)
samps <- sample.exact(knm, num.samps)
mrf.sample.plot(samps)

# Instantiate an empty model to fit:
psl <- make.crf(adj, n.states)
psl <- make.features(psl)
psl <- make.par(psl, 6)
psl$node.par[1,1,] <- 1
psl$node.par[2,1,] <- 2
psl$node.par[3,1,] <- 3

psl$edge.par[[1]][1,1,1] <- 4
psl$edge.par[[1]][2,2,1] <- 4
psl$edge.par[[2]][1,1,1] <- 5
psl$edge.par[[2]][2,2,1] <- 5
psl$edge.par[[3]][1,1,1] <- 6
psl$edge.par[[3]][2,2,1] <- 6
psl$edges

Delta.alpha.info <- delta.alpha(crf = psl, samples = samps, printQ = F)
Delta.alpha <- Delta.alpha.info$Delta.alpha

y <-c(samps[,1], samps[,2], samps[,3])
y
y[which(y==2)] <- 0
y

Ma <- glm(y ~ Delta.alpha[,1] + Delta.alpha[,2] + Delta.alpha[,3] +
              Delta.alpha[,4] + Delta.alpha[,5] + Delta.alpha[,6] - 1,
              family=binomial(link="logit"))
summary(Ma)

coef(Ma)
knm$par

psl$par <- as.numeric(coef(Ma))

out.pot2 <- make.pots(parms = psl$par,  crf = psl,  rescaleQ = T, replaceQ = T)
psl$node.pot
psl$edge.pot

# Node and edge beliefs:
psl.bel <- infer.exact(psl)
knm.bel <- infer.exact(knm)

psl.bel$node.bel
knm.bel$node.bel

psl.bel$edge.bel[[1]]
knm.bel$edge.bel[[1]]
psl.bel$edge.bel[[2]]
knm.bel$edge.bel[[2]]
psl.bel$edge.bel[[3]]
knm.bel$edge.bel[[3]]

# Configuration probabilities:
pot.info <- make.gRbase.potentials(psl, node.names = gp@nodes, state.nmes = c("1","2"))
pot.info

gR.dist.info    <- distribution.from.potentials(pot.info$node.potentials, pot.info$edge.potentials)
logZ            <- gR.dist.info$logZ
joint.dist.info <- as.data.frame(as.table(gR.dist.info$state.probs))
joint.dist.info

pot.info.knm <- make.gRbase.potentials(knm, node.names = gp@nodes, state.nmes = c("1","2"))
#pot.info.knm

gR.dist.info.knm    <- distribution.from.potentials(pot.info.knm$node.potentials, pot.info.knm$edge.potentials)
logZ.knm            <- gR.dist.info.knm$logZ
joint.dist.info.knm <- as.data.frame(as.table(gR.dist.info.knm$state.probs))
joint.dist.info.knm

psl.cp <- round(joint.dist.info[,4]*100, 1)
knm.cp <- round(joint.dist.info.knm[,4]*100, 1)
cbind(joint.dist.info[,c(2,3,1)], psl.cp, knm.cp)

sum(psl.cp)
sum(knm.cp)
