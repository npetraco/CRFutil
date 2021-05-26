library(CRFutil)
library(rstanarm)
library(rethinking)

grf.eq                <- ~A + B + A:B
AB                    <- make.empty.field(graph.eq = grf.eq, parameterization.typ = "standard", plotQ = F)
exp.theta.true        <- 1 # theta here is exp(omega)
# ******* theta true 11,22 is theta-times more likely than 12 or 21 ******** ?????
AB$node.pot[1,1] <- 3
AB$node.pot[2,1] <- 1
AB$edge.pot[[1]][1,1] <- exp.theta.true
AB$edge.pot[[1]][2,2] <- exp.theta.true

num.samps <- 500
samps     <- sample.exact(AB, num.samps)
dev.off()
mrf.sample.plot(samps)

marg.eg.info     <- marginal.edge.bayes.loglin(samps)
post.pots        <- marg.eg.info$rescaled.posterior.pots
marg.eg.bel.info <- marginal.edge.bels.bayes(post.pots)

marginal.edge.bayes.bels.plot(posterior.edge.belief.info = marg.eg.bel.info, type="X1|X2")
marginal.edge.bayes.bels.plot(posterior.edge.belief.info = marg.eg.bel.info, type="X2|X1")
