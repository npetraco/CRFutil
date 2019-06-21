library(CRFutil)
library(rstanarm)
library(rethinking)


# Using the same node potentials (0.8,0.2), (0.5,0.5), what does the re-scaled edge potential look like under independence??
num.samps <- 100

# Samples are independent
samps <- cbind(
  1+rbinom(n = num.samps, size = 1, prob = 0.2),
  1+rbinom(n = num.samps, size = 1, prob = 0.37)
)

dev.off()
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
head(junk$`Bel(X1,X2)`)
head(junk$`Bel(X1)`)
head(junk$`Bel(X2)`)
head(junk$`Bel(X1|X2)`)
head(junk$`Bel(X2|X1)`)
marginal.edge.bayes.bels.plot(junk, type="X1|X2", ymax=2500, edge.empirical.prob.info=NULL)

head(junk$`Bel(X1|X2)`)
head(junk$`Bel(X1)`)

names(junk)

