library(igraph)
library(gRbase)
library(CRFutil)

# Make up a random graph
g <- erdos.renyi.game(10, 0.4, typ="gnp")
dev.off()
plot(g)

# Get its adjacency matrix and genrate an MRF sample
adj <- as.matrix(as_adj(g))

f0 <- function(y){ as.numeric(c((y==1),(y==2)))}
rmod <- make.empty.field(adj.mat = adj, parameterization.typ = "standard", plotQ = T)
rmod$par <- runif(rmod$n.par,-1.5,1.5)
rmod$par # "true" theta
out.pot <- make.pots(parms = rmod$par,  crf = rmod,  rescaleQ = T, replaceQ = T)
rmod$edges    # Check against adjacency matrix
rmod$node.pot
rmod$edge.pot

# So now sample from the model as if we obtained an experimental sample:
num.samps <- 500
samps <- sample.exact(rmod, num.samps)
colnames(samps) <- 1:ncol(samps)
mrf.sample.plot(samps)


# Fit params in multiple ways:
# First get formula for graph
gf <- adj2formula(adj)
dev.off()
plot(ug(gf))

# Empirical
emp.dist <- fit_empirical(samps)
head(emp.dist)
plot(emp.dist[,11], typ="h")

# CRF MLE


# logistic regression (glm)
# Bayes logistic regression (Stan, loo, WAIC)
# log linear (loglin, glm)
# Bayes log linear (Poisson, Stan, loo, WAIC)
# Bayes zero-inflated (Stan, MODEL MATRIX?????)
# Bayes neg-binomial
# MLE zero-inflated, geg-binomial??

