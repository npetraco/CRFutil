library(CRFutil)
library(Rgraphviz)
library(MASS)
library(gRim)
library(rstanarm)
library(rethinking)

grf.eq                <- ~A + B + A:B
AB                    <- make.empty.field(graph.eq = grf.eq, parameterization.typ = "standard", plotQ = F)
theta.true            <- 2 # theta here is exp(omega)
# ******* theta true 11,22 is theta-times more likely than 12 or 21 ******** ?????
AB$edge.pot[[1]][1,1] <- theta.true
AB$edge.pot[[1]][2,2] <- theta.true

num.samps <- 500
samps     <- sample.exact(AB, num.samps)
mrf.sample.plot(samps)
head(samps)


megb <- marginal.edge.bayes.loglin(samps)
megb$coefficients
megb$model
megb$ses
megb$y
megb$stanfit
summary(megb)
pejunk <- as.matrix(megb)
hist(pejunk[,4], xlab="omega")

megb

pepotj <- exp(pejunk[,4])/exp(-pejunk[,4])
hist(pepotj)
abline(v=1)
median(pepotj)
sd(pepotj)

theta.prob.int <- HPDI(samples = pepotj, prob = 0.95)
names(theta.prob.int)
theta.prob.int
4 %in% theta.prob.int
PI(samples = pepotj, prob = 0.95)
