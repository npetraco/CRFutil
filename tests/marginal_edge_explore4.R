library(CRFutil)
library(MASS)

# Using the same node potentials (0.8,0.2), (0.5,0.5), what does the re-scaled edge potential look like under independence??
num.samps <- 10000

# Must label sample states 1 and 2
samps <- cbind(
  1+rbinom(n = num.samps, size = 1, prob = 0.2),
  1+rbinom(n = num.samps, size = 1, prob = 0.5)
)


mrf.sample.plot(samps)
#samps
junk.e <- marginal.edge.emp.pr(samps, printQ = T)

junk <- NULL
junk <- marginal.edge.mrf(samps)
dump.crf(crf = junk)
junk.b <- marginal.edge.bels(edge.mrf.obj = junk, node.names = c("X1","X2"), printQ = T)

rbind(junk.e$`Pr(X1)`, junk.e$`Pr(X1|X2)`, junk.b$`Bel(X1)`, junk.b$`Bel(X1|X2)`)
rbind(junk.e$`Pr(X2)`, junk.e$`Pr(X2|X1)`, junk.b$`Bel(X2)`, junk.b$`Bel(X2|X1)`)
make.pots(parms = junk$par, crf = junk, rescaleQ = T, replaceQ = F)[[2]][[1]]
exp(junk$par)[3]

chisq.test(junk.e$edge.contingency.tbl) # H0 X1 _||_ X2

# Bayes for CIs on edge param
marginal.edge.loglin(samps)

# ????
colnames(samps) <- c("X.1","X.2")

# Fold samples into contingency table for a look
X.cont <- xtabs(~., data=data.frame(samps))
X.cont

# Triangle model
colnames(samps)
grphf <- ~1:2
# Check the graph:
gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)


# Fit the Hojgaard model with loglm and gRim
fact.grphf <- ~X.1 + X.2 + X.1:X.2
loglin.mle <- loglm(fact.grphf, data=X.cont); # loglm from MASS which uses loglin in base
loglin.mle
summary(loglin.mle)

X.loglm.coefs <- coef(loglin.mle)
X.loglm.coefs

# Fit contingency table:
X.cont.fitted <- fitted(loglin.mle)
X.cont.fitted
X.cont

# Flatten fit contingency table into matrix-table form
X.counts.fitted <- data.frame(ftable(X.cont.fitted))
X.counts.fitted
sum(X.counts.fitted[,3]) # Equal to the sample size?

X.cont
