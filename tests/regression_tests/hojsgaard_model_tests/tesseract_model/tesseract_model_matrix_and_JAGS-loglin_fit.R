library(CRFutil)
library(CRFutilRcppComponents)
library(R2jags)


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
X.all <- sapply(1:ncol(X.all), function(xx){as.numeric(X.all[,xx])}) # Needed for C routines
freqs <- ftab[,ncol(ftab)]
head(X.all)
dim(X.all)
length(freqs)

plot(1:nrow(X.all), freqs,typ="h", ylab="config freqs",xlab="config #")



# Model Matrix with respect to graph
# Model Matrix with respect to graph
theta.pars     <- fix_node_and_edge_par(node_par = tess$node.par, edge_par = tess$edge.par) # Needed for C routines
tess$node.par.c <- theta.pars$node_par
tess$edge.par.c <- theta.pars$edge_par
M              <- compute_model_matrix(
  configs       = X.all,
  edge_mat      = tess$edges,
  node_par      = tess$node.par.c,
  edge_par      = tess$edge.par.c,
  num_params_in = tess$n.par)
dim(M)


dat <- list(
  p = ncol(M),
  N = nrow(M),
  y = freqs,
  Mmodl = M
)

fpth <- "inst/poisson_model.bug"
model_parameters <- c("alpha","theta")

# Run the model
# model_run <- jags(data = dat, #inits = init,
#                   parameters.to.save = model_parameters,
#                   model.file = fpth,
#                   n.chains = 1,
#                   n.iter = 10,
#                   n.burnin = 5,
#                   n.thin = 1)
#
# model_run <- jags.parallel(
#   data = dat, #inits = init,
#   parameters.to.save = model_parameters,
#   model.file = fpth,
#   n.chains = 4,
#   n.iter = 10,
#   n.burnin = 5,
#   n.thin = 1)

(55000 - 5000)/50 # got about 8000 samples with Stan in 2 hours

model_run <- jags.parallel(data = dat, #inits = init,
                 parameters.to.save = model_parameters,
                 model.file = fpth,
                 n.chains = 8,
                 n.iter = 55000,
                 n.burnin = 5000,
                 n.thin = 50)


model_run
