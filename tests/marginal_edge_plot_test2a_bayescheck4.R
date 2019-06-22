library(CRFutil)
library(rstanarm)
library(rethinking)

grf.eq <- ~X + Y + Z + Y:X + Y:Z + X:Z
XY <- make.empty.field(graph.eq = grf.eq, parameterization.typ = "standard", plotQ = T)

XY$node.pot[1,] <- c(0.8,0.2)       # X
XY$node.pot[2,] <- c(1-0.37,0.37)   # Y
XY$node.pot[3,] <- c(0.75,1 - 0.75) # Z

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

# eg.bels <- marginal.edge.bels.bayes(eg.post.pots)
# marginal.edge.bayes.bels.plot(eg.bels, type="X1|X2", ymax=2500) # Broken?????
# hist(eg.bels$`Bel(X2|X1)`[,1])
# head(eg.bels$`Bel(X2|X1)`)
#
# head(eg.bels$`Bel(X1,X2)`)
