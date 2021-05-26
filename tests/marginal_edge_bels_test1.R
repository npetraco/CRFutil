library(CRFutil)
library(rstanarm)
library(rethinking)

grf.eq <- ~X + Y + Z + X:Y + X:Z + Y:Z
XY <- make.empty.field(graph.eq = grf.eq, parameterization.typ = "standard", plotQ = T)
XY$edges

pX2 <- 0.2
pY2 <- 0.37
pZ2 <- 0.25
XY$node.pot[1,] <- c(1-pX2, pX2) # X
XY$node.pot[2,] <- c(1-pY2, pY2) # Y
XY$node.pot[3,] <- c(1-pZ2, pZ2) # Z

XY$edges
nn <- c("X","Y","Z")
XY$edge.pot[[1]][1,1] <- XY$edge.pot[[1]][2,2] <- 1.2; nn[XY$edges[1,] ]
XY$edge.pot[[2]][1,1] <- XY$edge.pot[[2]][2,2] <- 3;   nn[XY$edges[2,] ]
XY$edge.pot[[3]][1,1] <- XY$edge.pot[[3]][2,2] <- 0.8; nn[XY$edges[3,] ]

XY$node.pot
XY$edge.pot

num.samps <- 500
samps.all <- sample.exact(XY, num.samps)
mrf.sample.plot(samps.all)
head(samps.all)
samps <- samps.all[,XY$edges[1,]]; nn[XY$edges[1,] ] #
samps <- samps.all[,XY$edges[2,]]; nn[XY$edges[2,] ] #
samps <- samps.all[,XY$edges[3,]]; nn[XY$edges[3,] ] #


XY.emp.prs <- marginal.edge.emp.pr(samps) # Check sample margins. Marginals should be node pots and show independence
#XY.emp.prs
XY.emp.prs$`Pr(X1|X2)`
XY.emp.prs$`Pr(X1)`
XY.emp.prs$`Pr(X2|X1)`
XY.emp.prs$`Pr(X2)`

# Now fit a Bayes loglin model instead
marg.eg.info <- marginal.edge.bayes.loglin(samps)
eg.post.pots <- marg.eg.info$rescaled.posterior.pots

marginal.edge.bayes.plot(eg.post.pots, type = "logpot")
marginal.edge.bayes.plot(eg.post.pots, type = "pot")

eg.bels <- marginal.edge.bels.bayes.BROKEN(eg.post.pots) # Broken?????
hist(eg.bels$`Bel(X1|X2)`[,1], main=colnames(eg.bels$`Bel(X1|X2)`)[1])
hist(eg.bels$`Bel(X1|X2)`[,2], main=colnames(eg.bels$`Bel(X1|X2)`)[2])
hist(eg.bels$`Bel(X1|X2)`[,3], main=colnames(eg.bels$`Bel(X1|X2)`)[3])
hist(eg.bels$`Bel(X1|X2)`[,4], main=colnames(eg.bels$`Bel(X1|X2)`)[4])
XY.emp.prs$`Pr(X1|X2)`

hist(eg.bels$`Bel(X2|X1)`[,1], main=colnames(eg.bels$`Bel(X2|X1)`)[1])  # Broken?????
hist(eg.bels$`Bel(X2|X1)`[,2], main=colnames(eg.bels$`Bel(X2|X1)`)[2])
hist(eg.bels$`Bel(X2|X1)`[,3], main=colnames(eg.bels$`Bel(X2|X1)`)[3])
hist(eg.bels$`Bel(X2|X1)`[,4], main=colnames(eg.bels$`Bel(X2|X1)`)[4])
XY.emp.prs$`Pr(X2|X1)`

hist(eg.bels$`Bel(X1)`[,1], main=colnames(eg.bels$`Bel(X1)`)[1])
hist(eg.bels$`Bel(X1)`[,2], main=colnames(eg.bels$`Bel(X1)`)[2])
XY.emp.prs$`Pr(X1)`

eg.bels$`Bel(X2|X1)`[,1]
head(eg.bels$`Bel(X1,X2)`) # Problem....
head(eg.post.pots)


grf.eq.eg <- ~A + B + A:B
teg <- make.empty.field(graph.eq = grf.eq.eg, parameterization.typ = "standard", plotQ = F)
epn <- 1
eg.post.pots[epn,]
teg$node.pot[1,1] <- eg.post.pots[epn,][1]
teg$node.pot[2,1] <- eg.post.pots[epn,][2]
teg$node.pot

teg$edge.pot[[1]][1,1] <-  eg.post.pots[epn,][3]
teg$edge.pot[[1]][2,2] <-  eg.post.pots[epn,][3]
teg$edge.pot

dir.bels   <- infer.exact(teg)
indir.bels <- marginal.edge.bels(teg, printQ = F)

dir.bels$edge.bel[[1]]
indir.bels$`Bel(X1,X2)`
eg.bels$`Bel(X1,X2)`[epn,] # Problem...

flatten.marginal.edge.beliefs(indir.bels)
marginal.edge.bels.bayes.BROKEN(eg.post.pots[1:2,])

eg.post.pots[epn,]
> eg.post.pots[epn,]
pot.tau1  pot.tau2 pot.omega
10.125812  1.158845  1.319985
teg$node.pot
> teg$node.pot
[,1] [,2]
[1,] 10.125812    1
[2,]  1.158845    1
teg$edge.pot
> teg$edge.pot
[[1]]
[,1]     [,2]
[1,] 1.319985 1.000000
[2,] 1.000000 1.319985

marginal.edge.bayes.bels.plot(eg.bels, type="X1|X2", ymax=2500)
marginal.edge.bayes.bels.plot(eg.bels, type="X2|X1", ymax=2500)
