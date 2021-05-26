library(igraph)
library(gRbase)
library(CRFutil)
#library(rstanarm)
library(rstan)
library(MASS)

# Make up a random graph
g <- erdos.renyi.game(10, 0.6, typ="gnp")
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
plot(emp.dist[,ncol(emp.dist)], typ="h", xlab="configuration state#", ylab="Emp. Freq.")

# True model from true params
tru.dist <- fit_true(rmod)
head(tru.dist)
reordr.idxs <- reorder_configs(emp.dist[,1:(ncol(emp.dist)-1)], tru.dist[,1:(ncol(emp.dist)-1)])
tru.dist    <- tru.dist[reordr.idxs,]
plot(tru.dist[,ncol(emp.dist)], typ="h", xlab="configuration state#", ylab="True Freq.")

# CRF MLE
mle.dist    <- fit_mle(gf,samps,infer.exact, mag.grad.tol = 1e-2)
reordr.idxs <- reorder_configs(emp.dist[,1:(ncol(emp.dist)-1)], mle.dist[,1:(ncol(emp.dist)-1)])
mle.dist    <- mle.dist[reordr.idxs,]
plot(mle.dist[,ncol(emp.dist)], typ="h", xlab="configuration state#", ylab="MLE. Freq.")

# logistic regression (glm) CAUTION SLOW!!!!!!!!!!!!!!
logis.dist    <- fit_logistic(gf, samps)
reordr.idxs   <- reorder_configs(emp.dist[,1:(ncol(emp.dist)-1)], logis.dist[,1:(ncol(emp.dist)-1)])
logis.dist    <- logis.dist[reordr.idxs,]
plot(logis.dist[,ncol(emp.dist)], typ="h", xlab="configuration state#", ylab="Logistic Freq.")

# Bayes logistic regression (Stan, loo, WAIC)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
blogis.dist <- fit_bayes_logistic(gf, samps)
reordr.idxs <- reorder_configs(emp.dist[,1:10], blogis.dist[,1:10])
blogis.dist <- blogis.dist[reordr.idxs,]
plot(blogis.dist[,11], typ="h", xlab="configuration state#", ylab="Logistic Freq.")


# log linear (loglin, glm)
loglin.dist <- fit_loglinear(gf, samps)
reordr.idxs <- reorder_configs(emp.dist[,1:(ncol(emp.dist)-1)], loglin.dist[,1:(ncol(emp.dist)-1)])
loglin.dist <- loglin.dist[reordr.idxs,]
plot(loglin.dist[,ncol(emp.dist)], typ="h", xlab="configuration state#", ylab="Log-Linear Freq.")

# Bayes log linear (Poisson, Stan, loo, WAIC)  Doesn't work with this model matrix !!!!!!! FIX******************
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
bloglin.dist <- fit_bayes_loglinear(gf, samps)

reordr.idxs  <- reorder_configs(emp.dist[,1:10], bloglin.dist[,1:10])
bloglin.dist <- bloglin.dist[reordr.idxs,]
plot(bloglin.dist[,11], typ="h", xlab="configuration state#", ylab="Bayes Loglin Freq.")


# Bayes log linear2 (Poisson, Stan, loo, WAIC)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
model.c       <- stanc(file = "inst/poisson_model.stan", model_name = 'model')
sm            <- stan_model(stanc_ret = model.c, verbose = T)
bloglin2.dist <- fit_bayes_loglinear2(gf, samps, stan.model = sm)

reordr.idxs   <- reorder_configs(emp.dist[,1:10], bloglin2.dist[,1:10])
bloglin2.dist <- bloglin2.dist[reordr.idxs,]
plot(bloglin2.dist[,11], typ="h", xlab="configuration state#", ylab="Bayes Loglin2 Freq.")

# Bayes zero-inflated (Stan, MODEL MATRIX?????)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
model.c       <- stanc(file = "inst/zero_inflated_poisson_take3.stan", model_name = 'model')
sm            <- stan_model(stanc_ret = model.c, verbose = T)


bzipp.dist <- fit_bayes_zip(gf, samps, stan.model = sm, iter = 2000, chains = 4)
library(shinystan)
launch_shinystan(bzipp.dist)

si <- summary(bzipp.dist)
head(si$summary)
jneff <- si$summary[,9]
jR <- si$summary[,10]
min(jneff)
jneff[1]
min(jneff[-which(is.nan(jneff) == T)])
max(jR[-which(is.nan(jR) == T)])


sparams <- extract(bzipp.dist, permuted = TRUE)
colnames(sparams)
sparams$y_new[1:10]
length(sparams$y_new)
dim(sparams$beta)
length(rmod$par)
jtheta <- apply(sparams$beta, 2, median)
jtheta
jbeta <- extract(bzipp.dist, "beta")[[1]]
dim(jbeta)

jbzipp.fit   <- make.empty.field(graph.eq = gf, parameterization.typ = "standard")
jbzipp.fit$par   <- apply(extract(bzipp.dist,"beta_theta")[[1]], 2, median)

jout.potsx       <- make.pots(parms = bzipp.fit$par, crf = bzipp.fit, rescaleQ = T, replaceQ = T)

jpotentials.info    <- make.gRbase.potentials(bzipp.fit, node.names = colnames(samples), state.nmes = c("1","2"))
jdistribution.info  <- distribution.from.potentials(potentials.info$node.potentials, potentials.info$edge.potentials)
jjoint.distribution <- as.data.frame(as.table(distribution.info$state.probs))

# # Re-order columns to increasing order
# freq.idx    <- ncol(joint.distribution)
# node.nums   <- colnames(joint.distribution)[-freq.idx]
# node.nums   <- unlist(strsplit(node.nums, split = "X"))
# node.nums   <- node.nums[-which(node.nums == "")]
# node.nums   <- as.numeric(node.nums)
# col.reorder <- order(node.nums)
# joint.distribution <- joint.distribution[,c(col.reorder, freq.idx)]




jconfigs.and.counts <- as.data.frame(ftable(data.frame(samps)))
head(jconfigs.and.counts)
jbzipp.fit   <- make.empty.field(graph.eq = gf, parameterization.typ = "standard")

jM  <- compute.model.matrix(
  configs   = jconfigs.and.counts[,-ncol(jconfigs.and.counts)],
  edges.mat = jbzipp.fit$edges,
  node.par  = jbzipp.fit$node.par,
  edge.par  = jbzipp.fit$edge.par,
  ff        = f0)
jfreq <- jconfigs.and.counts[,ncol(jconfigs.and.counts)]




length(which(jconfigs.and.counts[,11] == 0))/nrow(jconfigs.and.counts)
ppi <- extract(bzipp.dist,"theta")[[1]]
hist(as.numeric(ppi))

bzipp.dist2 <- fit_bayes_zip(gf, samps, stan.model = sm, iter = 2000, chains = 4)



# Bayes neg-binomial
# MLE zero-inflated, neg-binomial??

# Assess difference from the true distribution
hist(mle.dist[,11]      - tru.dist[,11])
hist(emp.dist[,11]      - tru.dist[,11])
hist(logis.dist[,11]    - tru.dist[,11])
hist(blogis.dist[,11]   - tru.dist[,11])
hist(loglin.dist[,11]   - tru.dist[,11])
hist(bloglin.dist[,11]  - tru.dist[,11])
hist(bloglin2.dist[,11] - tru.dist[,11])
