library(igraph)
library(gRbase)
library(CRFutil)
#library(rstanarm)
library(rstan)
library(shinystan)
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


# Get frequency-counts of the configuration states:
configs.and.counts <- as.data.frame(ftable(data.frame(samps)))
freq.idx           <- ncol(configs.and.counts)
loc.freqs          <- configs.and.counts[,freq.idx]
loc.configs        <- configs.and.counts[,-freq.idx]

# Instantiate an empty model
bzipp.fit   <- make.empty.field(graph.eq = gf, parameterization.typ = "standard")
loc.f0      <- function(yy){ as.numeric(c((yy==1),(yy==2)))}

# Construct MRF model matrix.
print("Building model matrix")
loc.M  <- compute.model.matrix(
  configs   = loc.configs,
  edges.mat = bzipp.fit$edges,
  node.par  = bzipp.fit$node.par,
  edge.par  = bzipp.fit$edge.par,
  ff        = loc.f0)
print("Done with model matrix. Sorry it's slow...")
dim(loc.M)
length(rmod$par)
head(loc.M)
# For ZIP add intercept column
loc.M <- cbind(rep(1,nrow(loc.M)), loc.M)

loc.dat <- list(
  p = ncol(loc.M),
  N = nrow(loc.M),
  y = loc.freqs,
  X = loc.M,
  s = rep(10,ncol(loc.M)),
  s_theta = rep(10,ncol(loc.M))
)
#print(loc.dat)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
loc.model.c <- stanc(file = "inst/zero_inflated_poisson_take3b.stan", model_name = 'model')
loc.sm      <- stan_model(stanc_ret = loc.model.c, verbose = T)

loc.bfit <- sampling(loc.sm,
                     data    = loc.dat,
                     control = list(max_treedepth=20),
                     iter    = 5000,
                     thin    = 5,
                     chains  = 4)
fsinfo <- summary(loc.bfit)
dim(fsinfo$summary)
neffs <- fsinfo$summary[,9]
rs    <- fsinfo$summary[,10]
hist(neffs)
min(neffs,na.rm = T)
max(rs,na.rm = T)

pars <- c("beta_lam","beta_ppi","lp__")
print(loc.bfit, pars = pars)

pars2 <- c("lambda","ppi")
print(loc.bfit, pars = pars2)

pairs(loc.bfit, pars = pars)

launch_shinystan(loc.bfit)
check_hmc_diagnostics(loc.bfit)

# First beta is the intercept
bzipp.fit$par   <- apply(extract(loc.bfit,"beta_lam")[[1]], 2, median)[-1]
bzipp.fit$par
rmod$par

out.potsx       <- make.pots(parms = bzipp.fit$par, crf = bzipp.fit, rescaleQ = F, replaceQ = T)

potentials.info    <- make.gRbase.potentials(bzipp.fit, node.names = colnames(samps), state.nmes = c("1","2"))
distribution.info  <- distribution.from.potentials(potentials.info$node.potentials, potentials.info$edge.potentials)
joint.distribution <- as.data.frame(as.table(distribution.info$state.probs))

# Re-order columns to increasing order
freq.idx    <- ncol(joint.distribution)
node.nums   <- colnames(joint.distribution)[-freq.idx]
node.nums   <- unlist(strsplit(node.nums, split = "X"))
node.nums   <- node.nums[-which(node.nums == "")]
node.nums   <- as.numeric(node.nums)
col.reorder <- order(node.nums)
jnt.dist  <- joint.distribution[,c(col.reorder, freq.idx)]

reordr.idxs <- reorder_configs(emp.dist[,1:(ncol(emp.dist)-1)], jnt.dist[,1:(ncol(emp.dist)-1)])
jnt.dist <- jnt.dist[reordr.idxs,]
plot(jnt.dist[,ncol(emp.dist)], typ="h", xlab="configuration state#", ylab="ZIP Freq.")

hist(jnt.dist[,ncol(emp.dist)] - tru.dist[,ncol(emp.dist)])
hist(emp.dist[,ncol(emp.dist)] - tru.dist[,ncol(emp.dist)])

plot(tru.dist[,ncol(emp.dist)], typ="h", xlab="configuration state#", ylab="True Freq.")
plot(emp.dist[,ncol(emp.dist)], typ="h", xlab="configuration state#", ylab="Emp. Freq.")
plot(jnt.dist[,ncol(emp.dist)], typ="h", xlab="configuration state#", ylab="ZIP Freq.")

#save(samps,bzipp.fit,gf,loc.bfit, jnt.dist, file = "tests/regression_tests/hojsgaard_model_tests/randomg_zip.RData")

loglam <- extract(loc.bfit, "lambda")[[1]]
lam    <- exp(loglam)
dim(lam)
lam.med <- apply(lam,2,median)
prc.med <- lam.med/nrow(loc.configs)
prc.med[reordr.idxs]

plot(emp.dist[,ncol(emp.dist)], typ="h", xlab="configuration state#", ylab="Emp. Freq.")
plot(jnt.dist[,ncol(emp.dist)], typ="h", xlab="configuration state#", ylab="ZIP Freq.")
plot(prc.med, typ="h", xlab="configuration state#", ylab="ZIP Freq.")
sum(prc.med)
sum(jnt.dist[,ncol(emp.dist)])

# ????????
logitppi     <- extract(loc.bfit, "ppi")[[1]]
ppi <- exp(logitppi)/(1+exp(logitppi))

ppi.med <- apply(ppi,2,median)
plot(ppi.med, typ="h", xlab="configuration state#", ylab="ZIP ppi")
#logit()
#intercept
