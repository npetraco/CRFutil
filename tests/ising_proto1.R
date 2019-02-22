library(CRFutil)
library(rstan)
library(shinystan)
library(coda)

# Model:
edgs <-  rbind(
  c(1,2),
  c(1,5),
  c(2,6),
  c(2,3),
  c(3,7),
  c(3,4),
  c(4,8),
  c(5,6),
  c(5,9),
  c(6,7),
  c(6,10),
  c(7,8),
  c(7,11),
  c(8,12),
  c(9,13),
  c(9,10),
  c(10,11),
  c(10,14),
  c(11,12),
  c(11,15),
  c(12,16),
  c(13,14),
  c(14,15),
  c(15,16)
)

adj <- edges2adj(edgs, plotQ = T)
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}
n.states <- 2


knm <- make.empty.field(
  graph.eq             = NULL,
  adj.mat              = adj,
  parameterization.typ = "ising1",
  node.par             = NULL,
  edge.par             = NULL,
  plotQ                = T)
dump.crf(knm)
knm$edges
knm$par
knm$n.par
knm$node.par
knm$node.pot

knm$par <- -1.145
out.pot <- make.pots(parms = knm$par,  crf = knm,  rescaleQ = T, replaceQ = T)
knm$edge.par
knm$edge.pot


# So now sample from the model as if we obtained an experimental sample:
num.samps <- 1000
#set.seed(1)
samps <- sample.exact(knm, num.samps)
mrf.sample.plot(samps)

psl <- make.empty.field(
  graph.eq             = NULL,
  adj.mat              = adj,
  parameterization.typ = "ising1",
  node.par             = NULL,
  edge.par             = NULL,
  plotQ                = F)

Delta.alpha.info <- delta.alpha(crf = psl, samples = samps, printQ = F)
Delta.alpha <- Delta.alpha.info$Delta.alpha
Delta.alpha

y <-as.numeric(samps)
y
y[which(y==2)] <- 0
y


Ma <- glm(y ~ Delta.alpha[,1] - 1, family=binomial(link="logit"))
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


# Config probs:
pot.info.knm <- make.gRbase.potentials(knm, node.names = 1:16, state.nmes = c("1","2"))
gR.dist.info.knm    <- distribution.from.potentials(pot.info.knm$node.potentials, pot.info.knm$edge.potentials)
logZ.knm            <- gR.dist.info.knm$logZ
joint.dist.info.knm <- as.data.frame(as.table(gR.dist.info.knm$state.probs))
dim(joint.dist.info.knm)
head(joint.dist.info.knm,10)
plot(joint.dist.info.knm[,17], typ="h")


pot.info.psl <- make.gRbase.potentials(psl, node.names = 1:16, state.nmes = c("1","2"))
gR.dist.info.psl    <- distribution.from.potentials(pot.info.psl$node.potentials, pot.info.psl$edge.potentials)
logZ.psl            <- gR.dist.info.psl$logZ
joint.dist.info.psl <- as.data.frame(as.table(gR.dist.info.psl$state.probs))
dim(joint.dist.info.psl)
head(joint.dist.info.psl,10)
plot(joint.dist.info.psl[,17], typ="h")


psl.cp <- round(joint.dist.info.psl[,17]*100, 6)    # MLE logistic
knm.cp <- round(joint.dist.info.knm[,17]*100, 6)    # True
head(cbind(knm.cp, psl.cp), 10)

sum(psl.cp)
sum(knm.cp)


samps
library(reshape2)
library(ggplot2)

#snum <- sample(1:nrow(samps),1)
snum <- which(joint.dist.info.knm[,17]>0.1)[1] # High bars in the joint show congigs that should come up alot, and hence have a high chance of being seen multiple times on random host objects
#x=t(matrix(samps[snum,], ncol=4))

x=t(matrix(as.numeric(joint.dist.info.knm[snum,1:16]), ncol=4))
x
x1=melt(x)
names(x1)=c("x","y","color")
x1$color=factor(x1$color>1)
levels(x1$color)=c("1","2")
qplot(x, y, fill=color, data=x1,geom='tile')
snum

x
which(joint.dist.info.knm[,17]>0.1)[2]
dim(samps)

16^2
2^256
# Demo for increasing size but still small fields
# 16x16 Can't get Prs for all these configurations, but could be get Z ??
2^64

# loop ove nodes (last node??)
# determine node type
# input edges for type along with their designation (horr, vert, diag)

