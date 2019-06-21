library(CRFutil)
library(rstanarm)
library(rethinking)


# Using the same node potentials (0.8,0.2), (0.5,0.5), what does the re-scaled edge potential look like under independence??
num.samps <- 5

samps <- cbind(
  1+rbinom(n = num.samps, size = 1, prob = 0.2),
  1+rbinom(n = num.samps, size = 1, prob = 0.5)
)

dev.off()
mrf.sample.plot(samps)

marg.eg.info     <- marginal.edge.bayes.loglin(samps)
post.pots        <- marg.eg.info$rescaled.posterior.pots
marg.eg.bel.info <- marginal.edge.bels.bayes(post.pots)

marginal.edge.bayes.bels.plot(posterior.edge.belief.info = marg.eg.bel.info, type="X1|X2")
marginal.edge.bayes.bels.plot(posterior.edge.belief.info = marg.eg.bel.info, type="X2|X1")

