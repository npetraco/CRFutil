library(CRFutil)
#library(Rgraphviz)
#library(MASS)
#library(gRim)
library(rstanarm)
library(rethinking)

grf.eq                <- ~A + B + A:B
AB                    <- make.empty.field(graph.eq = grf.eq, parameterization.typ = "standard", plotQ = F)
theta.true            <- 1.5 # theta here is exp(omega)
# ******* theta true 11,22 is theta-times more likely than 12 or 21 ******** ?????
AB$node.pot[1,1] <- 3
AB$edge.pot[[1]][1,1] <- theta.true
AB$edge.pot[[1]][2,2] <- theta.true

num.samps <- 500
samps     <- sample.exact(AB, num.samps)
mrf.sample.plot(samps)
head(samps)

megb <- marginal.edge.bayes.loglin(samps)
megb$raw.theta
megb$theta.par.meds
ppots <- megb$rescaled.posterior.pots
head(ppots)

colnames(ppots)
int.prob.lvl <- 0.95
pptau1.int <- HPDI(ppots[,1], prob = int.prob.lvl)
pptau2.int <- HPDI(ppots[,2], prob = int.prob.lvl)
ppomeg.int <- HPDI(ppots[,3], prob = int.prob.lvl)

hist(ppots[,1], xlab="exp(tau1)");abline(v=1, lwd=4, col="green")
points(pptau1.int, c(0,0), pch=16,col="blue")
pptau1.int

hist(ppots[,2], xlab="exp(tau2)");abline(v=1, lwd=3, col="green")
points(pptau2.int, c(0,0), pch=16,col="blue")
pptau2.int

hist(ppots[,3], xlab="exp(omega)");abline(v=1, lwd=3, col="green")
points(ppomeg.int, c(0,0), pch=16,col="blue")
ppomeg.int

junk <- marginal.edge.bels.bayes(ppots)

hist(junk$"Bel(X1)"[,1], xlab="Bel(X1=1)", main="")
hist(junk$"Bel(X1)"[,2], xlab="Bel(X1=2)", main="")

hist(junk$"Bel(X2)"[,1], xlab="Bel(X2=1)", main="")
hist(junk$"Bel(X2)"[,2], xlab="Bel(X2=2)", main="")

head(junk$"Bel(X1,X2)")
hist(junk$"Bel(X1,X2)"[,1], xlab="Bel(X1=1, X2=1)", main="")
hist(junk$"Bel(X1,X2)"[,2], xlab="Bel(X1=2, X2=1)", main="")
hist(junk$"Bel(X1,X2)"[,3], xlab="Bel(X1=1, X2=2)", main="")
hist(junk$"Bel(X1,X2)"[,4], xlab="Bel(X1=2, X2=2)", main="")


head(junk$"Bel(X1|X2)")
hist(junk$"Bel(X1|X2)"[,1], xlab="Bel(X1=1 | X2=1)", main="", xlim=c(0,1), ylim=c(0,2500))
abline(v=median(junk$"Bel(X1)"[,"X1=1"]), lwd=4)
median(junk$"Bel(X1)"[,"X1=1"])
par(new=T)
hist(junk$"Bel(X1)"[,"X1=1"], main="Bel(X1=1)", xlab="", xlim=c(0,1), ylim=c(0,2500))


hist(junk$"Bel(X1|X2)"[,3], xlab="Bel(X1=1 | X2=2)", main="", xlim=c(0,1), ylim=c(0,2500))
     abline(v=median(junk$"Bel(X1)"[,"X1=1"]), lwd=4)
              median(junk$"Bel(X1)"[,"X1=1"])
par(new=T)
hist(junk$"Bel(X1)"[,"X1=1"], main="Bel(X1=1)", xlab="", xlim=c(0,1), ylim=c(0,2500))

hist(junk$"Bel(X1|X2)"[,2], xlab="Bel(X1=2 | X2=1)", main="", xlim=c(0,1), ylim=c(0,2500))
     abline(v=median(junk$"Bel(X1)"[,"X1=2"]), lwd=4)
              median(junk$"Bel(X1)"[,"X1=2"])
par(new=T)
hist(junk$"Bel(X1)"[,"X1=2"], main="Bel(X1=1)", xlab="", xlim=c(0,1), ylim=c(0,2500))

hist(junk$"Bel(X1|X2)"[,4], xlab="Bel(X1=2 | X2=2)", main="", xlim=c(0,1), ylim=c(0,2500))
abline(v=median(junk$"Bel(X1)"[,"X1=2"]), lwd=4)
median(junk$"Bel(X1)"[,"X1=2"])
par(new=T)
hist(junk$"Bel(X1)"[,"X1=2"], main="Bel(X1=1)", xlab="", xlim=c(0,1), ylim=c(0,2500))

head(junk$"Bel(X1|X2)")
hist(junk$"Bel(X1|X2)"[,1], xlab="Bel(X1=1 | X2=1)", main="")
hist(junk$"Bel(X1|X2)"[,2], xlab="Bel(X1=2 | X2=1)", main="")
hist(junk$"Bel(X1|X2)"[,3], xlab="Bel(X1=1 | X2=2)", main="")
hist(junk$"Bel(X1|X2)"[,4], xlab="Bel(X1=2 | X2=2)", main="")

head(junk$"Bel(X2|X1)")
hist(junk$"Bel(X2|X1)"[,1], xlab="Bel(X2=1 | X1=1)", main="")
hist(junk$"Bel(X2|X1)"[,2], xlab="Bel(X2=2 | X1=1)", main="")
hist(junk$"Bel(X2|X1)"[,3], xlab="Bel(X2=1 | X1=2)", main="")
hist(junk$"Bel(X2|X1)"[,4], xlab="Bel(X2=2 | X1=2)", main="")

