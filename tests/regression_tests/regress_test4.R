library(CRFutil)
library(rstan)
library(shinystan)

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

set.seed(6)
knm$par <- runif(6,-1.5,1.1)
knm$par # "true" theta
out.pot <- make.pots(parms = knm$par,  crf = knm,  rescaleQ = T, replaceQ = T)

# So now sample from the model as if we obtained an experimental sample:
num.samps <- 25
set.seed(1)
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

#----------
psl2 <- make.crf(adj, n.states)
psl2 <- make.features(psl2)
psl2 <- make.par(psl2, 6)
psl2$node.par[1,1,] <- 1
psl2$node.par[2,1,] <- 2
psl2$node.par[3,1,] <- 3

psl2$edge.par[[1]][1,1,1] <- 4
psl2$edge.par[[1]][2,2,1] <- 4
psl2$edge.par[[2]][1,1,1] <- 5
psl2$edge.par[[2]][2,2,1] <- 5
psl2$edge.par[[3]][1,1,1] <- 6
psl2$edge.par[[3]][2,2,1] <- 6
#----------

Delta.alpha.info <- delta.alpha(crf = psl, samples = samps, printQ = F)
Delta.alpha <- Delta.alpha.info$Delta.alpha

y <-c(samps[,1], samps[,2], samps[,3])
y
y[which(y==2)] <- 0
y

dat <- list(
  N=length(y),
  X1=Delta.alpha[,1],
  X2=Delta.alpha[,2],
  X3=Delta.alpha[,3],
  X4=Delta.alpha[,4],
  X5=Delta.alpha[,5],
  X6=Delta.alpha[,6],
  y = y
)

model.c <- stanc(file = "/Users/npetraco/latex/papers/dust/steph_diss/CRFutil/tests/regression_tests/logistic_model1.stan", model_name = 'model')
sm <- stan_model(stanc_ret = model.c, verbose = T)
bfit <- sampling(sm, data = dat, iter=10000, thin = 1, chains = 4)

beta <- extract(bfit,"beta")[[1]]
dim(beta)
plot(beta[,1], typ="l", main="Trace")
hist(beta[,1], bre=80, probability = T)
acf(beta[,1])
mean(beta[,1])
median(beta[,1])

launch_shinystan(bfit)

Ma <- glm(y ~ Delta.alpha[,1] + Delta.alpha[,2] + Delta.alpha[,3] +
              Delta.alpha[,4] + Delta.alpha[,5] + Delta.alpha[,6] - 1,
              family=binomial(link="logit"))
summary(Ma)


c(mean(beta[,1]), mean(beta[,2]), mean(beta[,3]), mean(beta[,4]), mean(beta[,5]), mean(beta[,6]))
c(median(beta[,1]), median(beta[,2]), median(beta[,3]), median(beta[,4]), median(beta[,5]), median(beta[,6]))
knm$par

psl$par <- c(median(beta[,1]), median(beta[,2]), median(beta[,3]), median(beta[,4]), median(beta[,5]), median(beta[,6]))
psl2$par <- as.numeric(coef(Ma))

out.pot2 <- make.pots(parms = psl$par,  crf = psl,  rescaleQ = T, replaceQ = T)
psl$node.pot
psl$edge.pot

out.pot2a <- make.pots(parms = psl2$par,  crf = psl2,  rescaleQ = T, replaceQ = T)
psl2$node.pot
psl2$edge.pot


# Node and edge beliefs:
psl.bel  <- infer.exact(psl)
psl2.bel <- infer.exact(psl2)
knm.bel  <- infer.exact(knm)

psl.bel$node.bel
psl2.bel$node.bel
knm.bel$node.bel

psl.bel$edge.bel[[1]]
psl2.bel$edge.bel[[1]]
knm.bel$edge.bel[[1]]


psl.bel$edge.bel[[2]]
psl2.bel$edge.bel[[2]]
knm.bel$edge.bel[[2]]


psl.bel$edge.bel[[3]]
psl2.bel$edge.bel[[3]]
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


pot.info2 <- make.gRbase.potentials(psl2, node.names = gp@nodes, state.nmes = c("1","2"))
pot.info2

gR.dist.info2    <- distribution.from.potentials(pot.info2$node.potentials, pot.info2$edge.potentials)
logZ2            <- gR.dist.info2$logZ
joint.dist.info2 <- as.data.frame(as.table(gR.dist.info2$state.probs))
joint.dist.info2

pot.info.knm <- make.gRbase.potentials(knm, node.names = gp@nodes, state.nmes = c("1","2"))


gR.dist.info.knm    <- distribution.from.potentials(pot.info.knm$node.potentials, pot.info.knm$edge.potentials)
logZ.knm            <- gR.dist.info.knm$logZ
joint.dist.info.knm <- as.data.frame(as.table(gR.dist.info.knm$state.probs))
joint.dist.info.knm

psl.cp  <- round(joint.dist.info[,4]*100, 1)     # Bayes logistic
psl2.cp <- round(joint.dist.info2[,4]*100, 1)    # MLE logistic
knm.cp  <- round(joint.dist.info.knm[,4]*100, 1) # True
cbind(joint.dist.info[,c(2,3,1)], psl.cp, psl2.cp, knm.cp)

sum(psl.cp)
sum(knm.cp)
