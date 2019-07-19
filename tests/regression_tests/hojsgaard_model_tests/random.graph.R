library(igraph)
library(gRbase)
library(CRFutil)

# Make up a random graph
g <- erdos.renyi.game(10, 0.4, typ="gnp")
dev.off()
plot(g)

# Get its adjacency matrix and genrate an MRF sample
adj <- as.matrix(as_adj(g))

f0       <- function(y){ as.numeric(c((y==1),(y==2)))}
rmod     <- make.empty.field(adj.mat = adj, parameterization.typ = "standard", plotQ = F)

# "true" theta
rmod$par <- runif(rmod$n.par,-1.5,1.5)
rmod$par

# Make true pots from true theta
out.pot <- make.pots(parms = rmod$par,  crf = rmod,  rescaleQ = T, replaceQ = T)
rmod$edges
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
dev.off()
plot(emp.dist[,11], typ="h", xlab="configuration state#", ylab="Emp. Freq.")

# True model from true params
tru.dist <- fit_true(rmod)
head(tru.dist)
reordr.idxs <- reorder_configs(emp.dist[,1:10], tru.dist[,1:10])
tru.dist    <- tru.dist[reordr.idxs,]
plot(tru.dist[,11], typ="h", xlab="configuration state#", ylab="True Freq.")

# CRF MLE
mle.dist    <- fit_mle(gf,samps,infer.exact)
reordr.idxs <- reorder_configs(emp.dist[,1:10], mle.dist[,1:10])
mle.dist    <- mle.dist[reordr.idxs,]
plot(mle.dist[,11], typ="h", xlab="configuration state#", ylab="MLE. Freq.")

# logistic regression (glm)
logis.dist    <-fit_logistic(gf, samps)
reordr.idxs   <- reorder_configs(emp.dist[,1:10], logis.dist[,1:10])
logis.dist    <- logis.dist[reordr.idxs,]
plot(logis.dist[,11], typ="h", xlab="configuration state#", ylab="Logistic Freq.")


# Bayes logistic regression (Stan, loo, WAIC)
# log linear (loglin, glm)
# Bayes log linear (Poisson, Stan, loo, WAIC)
# Bayes zero-inflated (Stan, MODEL MATRIX?????)
# Bayes neg-binomial
# MLE zero-inflated, geg-binomial??

# Assess difference from the true distribution
hist(mle.dist[,11] - tru.dist[,11])
hist(emp.dist[,11] - tru.dist[,11])
hist(logis.dist[,11] - tru.dist[,11])
