# loglin gRim, model matrix loglin, bayes logistic, mle logistic, CRF fit

# loglin gRim
library(gRim)
library(MASS)

# Triangle model
grphf <- ~1:2+1:3+2:3
# Check the graph:
gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)


#         FIX
# First fit the Hojgaard model with loglm and gRim
ll2 <- loglm(grphf, data=X); # loglm from MASS which uses loglin in base
ll2
X.loglm.coefs <- coef(ll2)
X.loglm.coefs

# Look at emp relative freqs vs the fitted relative freqs:
X.fitted <- fitted(ll2)
X.Prob.fitted <- X.fitted/sum(X.fitted)
ftable(X.Prob.fitted)
data.frame(as.data.frame(as.table(X.Prob)), as.data.frame(as.table(X.Prob.fitted)))
