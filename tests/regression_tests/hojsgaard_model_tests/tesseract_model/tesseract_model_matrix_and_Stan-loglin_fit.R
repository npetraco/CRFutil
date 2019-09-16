library(rstan)
library(shinystan)
library(CRFutil)


# Tesseract field model:
grphf <-
 ~X.1:X.2   + X.1:X.4   + X.1:X.5  + X.1:X.13 +
  X.2:X.3   + X.2:X.6   + X.2:X.14 +
  X.3:X.4   + X.3:X.7   + X.3:X.15 +
  X.4:X.8   + X.4:X.16  +
  X.5:X.6   + X.5:X.8   + X.5:X.9 +
  X.6:X.7   + X.6:X.10  +
  X.7:X.8   + X.7:X.11  +
  X.8:X.12  +
  X.9:X.10  + X.9:X.12  + X.9:X.13 +
  X.10:X.11 + X.10:X.14 +
  X.11:X.12 + X.11:X.15 +
  X.12:X.16 +
  X.13:X.14 + X.13:X.16 +
  X.14:X.15 +
  X.15:X.16

adj <- ug(grphf, result="matrix")
adj
node2nme <- data.frame(1:ncol(adj),colnames(adj))
colnames(node2nme) <- c("node","name")
node2nme

gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)
plot(gp)

# Make up random potentials/sample and return a CRF-object
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}
tess <- make.empty.field(adj.mat = adj, parameterization.typ = "standard", plotQ = T)
set.seed(87)
tess$par <- runif(tess$n.par,-1.5,1.5)
tess$par # "true" theta
out.pot <- make.pots(parms = tess$par,  crf = tess,  rescaleQ = F, replaceQ = T)
tess$edges
tess$node.pot
tess$edge.pot

# So now sample from the model as if we obtained an experimental sample:
num.samps <- 500
set.seed(1)
samps <- sample.exact(tess, num.samps)
mrf.sample.plot(samps)

node.names      <- colnames(adj)
node.names
colnames(samps) <- node.names
head(samps)

# Make state count freq table from samples and compute frequencies of all possible state configs.
# State configs not observed will have 0 freq
ftab <- data.frame(ftable(data.frame(samps)))
X.all <- ftab[,1:ncol(samps)]
freqs <- ftab[,ncol(ftab)]

X.all
class(X.all)
dim(X.all)
freqs
hist(freqs)
sum(freqs)


# Model Matrix with respect to graph
M <- compute.model.matrix(configs=X.all, edges.mat=tess$edges, node.par = tess$node.par, edge.par = tess$edge.par, ff = f0)
#M
dim(M)

# M.red <- compute.model.matrix(configs=X.red, edges.mat=tess$edges, node.par = tess$node.par, edge.par = tess$edge.par, ff = f0)
# dim(M.red)
# M.red

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Fit log(lam) = alp + M theta
fpth <- "C:/Users/aliso/codes/CRFutil/tests/regression_tests/hojsgaard_model_tests/triangle_model/"
#fpth <- "/home/npetraco/codes/R/CRFutil/tests/regression_tests/hojsgaard_model_tests/triangle_model/"
model.c <- stanc(file = paste0(fpth,"vanalla.poisson.regression2b.stan"), model_name = 'model')
sm <- stan_model(stanc_ret = model.c, verbose = T)


dat <- list(
  p = ncol(M),
  N = nrow(M),
  y = freqs,
  Mmodl = M
)
dat

bfit <- sampling(sm, data = dat, iter=2000, thin = 1, chains = 8)
save(file = "bfit.RData",bfit)
save(file = "M.RData",M)
save(file = "samps.RData",samps)

# Get the parameter estimtes:
theta <- extract(bfit, "theta")[[1]]
alpha <- extract(bfit, "alpha")[[1]]


# Quick checks:
th <- apply(theta, MARGIN = 2, FUN = median)
length(th)
length(tess$par)
dim(M)
Emat <- M%*%th
dim(Emat)
alpv <- median(alpha)*rep(1,nrow(M))
length(alpv)
logl <- alpv+Emat
plot(exp(logl),typ="h")
prs <- exp(logl)/nrow(samps)
plot(prs,typ="h")
sum(prs)

# Extract the graph's partition function:
#N    <- nrow(samps)
N <- nrow(samps)
logZ <- log(N) - alpha
hist(logZ)



# Now using the (medians) posterior theta, estimate partition
# function with Belief propegation
llm2     <- make.empty.field(adj.mat = adj, parameterization.typ = "standard")
llm2$par <- apply(theta, MARGIN = 2, FUN = median)
out.pot  <- make.pots(parms = llm2$par,  crf = llm2,  rescaleQ = F, replaceQ = T)

pot.info.llm2.model     <- make.gRbase.potentials(llm2, node.names = c("X.1", "X.2", "X.3","X.4",
                                                                       "X.5", "X.6", "X.7","X.8",
                                                                       "X.9", "X.10", "X.11","X.12",
                                                                       "X.13", "X.14", "X.15","X.16"), state.nmes = c("1","2"))
gR.dist.info.llm2.model <- distribution.from.potentials(pot.info.llm2.model$node.potentials, pot.info.llm2.model$edge.potentials)
logZ.llm2.model         <- gR.dist.info.llm2.model$logZ

logZ.llm2.model # log partition function from BP with a theta (from regression)
mean(logZ)      # log partition function direct from regression
median(logZ)

joint.dist.info.llm2.model <- as.data.frame(as.table(gR.dist.info.llm2.model$state.probs))
dim(joint.dist.info.llm2.model)
prs.bp <- joint.dist.info.llm2.model[,17]
sum(prs.bp)
plot(prs.bp,typ="h")
# Order configs!!!!


joint.dist.info.llm2.model <- joint.dist.info.llm2.model[,c(2,3,1,4)]
joint.dist.info.llm2.model[,4] <- joint.dist.info.llm2.model[,4] * 100
joint.dist.info.llm2.model


# OLD BELOW:

# Can't handle sampling on whole model
# Just try non-zero freques
# nz.idxs <- which(freqs>0)
# X.red <- X.all[nz.idxs,]
# freqs.red <- freqs[nz.idxs]
# dim(X.red)
# length(freqs.red)
# hist(freqs.red)
# freqs.red
# max(freqs.red)
# X.red[which(freqs.red == max(freqs.red)), ]


fit.smy <- summary(bfit)$summary
rhats <- fit.smy[,"Rhat"]
neffs <- fit.smy[,"n_eff"]
min(neffs)
max(rhats)

launch_shinystan(bfit)

theta <- extract(bfit, "theta")[[1]]
alpha <- extract(bfit, "alpha")[[1]]

hist(alpha)
NN <- nrow(samps)
logZ <- log(NN) - alpha
hist(logZ)

lambda <- extract(bfit, "lambda")[[1]]
config.probs <- colMeans(lambda)/NN      # Approx config probs
sum(config.probs) # Add about to 1?

max(config.probs)

# Above assumes configs dropped (zero freqs) had 0 prob!!
th <- colMeans(theta)
th
tess$par

# Check tess config probs vs lambda derived config probs
pot.info.tess.model        <- make.gRbase.potentials(tess, node.names = node.names, state.nmes = c("1","2"))
gR.dist.info.tess.model    <- distribution.from.potentials(pot.info.tess.model$node.potentials, pot.info.tess.model$edge.potentials)

gR.dist.info.tess.model$state.probs
gR.dist.info.tess.model$logZ


joint.dist.info.tess.model <- as.data.frame(as.table(gR.dist.info.tess.model$state.probs))
head(joint.dist.info.tess.model)
dim(joint.dist.info.tess.model)


col.reorder.idxs <- sapply(1:16, function(xx){which(colnames(joint.dist.info.tess.model)[-17] == node.names[xx])})


joint.dist.info.tess.model <- joint.dist.info.tess.model[,c(col.reorder.idxs,17)]
head(joint.dist.info.tess.model)
joint.dist.info.tess.model[,17] <- joint.dist.info.tess.model[,17] * 100


sample.config.idxs <- sapply(1:nrow(X.red), function(xx){row.match(X.red[xx,], joint.dist.info.tess.model[,1:16])})

data.frame(joint.dist.info.tess.model[sample.config.idxs,], config.probs)

sum(joint.dist.info.tess.model[sample.config.idxs,17]) # WOW! We miss alot by cutting out the configs not observed in the sample


