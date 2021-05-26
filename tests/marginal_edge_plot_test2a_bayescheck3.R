library(CRFutil)
library(rstanarm)
library(rethinking)

grf.eq <- ~X + Y + X:Y
XY <- make.empty.field(graph.eq = grf.eq, parameterization.typ = "standard", plotQ = F)

XY$node.pot[1,] <- c(0.8,0.2)
XY$node.pot[2,] <- c(1-0.37,0.37)
XY$edge.pot[[1]][1,1] <- XY$edge.pot[[1]][2,2] <- 1 # ******* theta true 11,22 is theta-times more likely than 12 or 21
XY$node.pot
XY$edge.pot

num.samps <- 500
samps     <- sample.exact(XY, num.samps)
mrf.sample.plot(samps)
XY.emp.prs <- marginal.edge.emp.pr(samps) # Check sample margins. Marginals should be node pots and show independence
XY.emp.prs

# Now fit a Bayes loglin model instead
marg.eg.info <- marginal.edge.bayes.loglin(samps)
eg.post.pots <- marg.eg.info$rescaled.posterior.pots

# Look at posterior of exp(omega) = pot.omega. It should comfortably cover 1 (0 for log) since
# the sample is independent
hist(eg.post.pots[,3])
hist(log(eg.post.pots[,3]))

# XXXXX
junk <- marginal.edge.bels.bayes(eg.post.pots)
#head(junk$`Bel(X1,X2)`)
head(junk$`Bel(X1|X2)`)
head(junk$`Bel(X1)`)
head(junk$`Bel(X2|X1)`)
head(junk$`Bel(X2)`)

marginal.edge.bayes.bels.plot(junk, type="X2|X1", ymax=3500, edge.empirical.prob.info=NULL)

head(eg.post.pots)
marginal.edge.bayes.plot(eg.post.pots, type = "logpot")
