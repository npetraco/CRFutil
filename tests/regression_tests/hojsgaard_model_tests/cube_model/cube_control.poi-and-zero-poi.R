library(rstan)
library(shinystan)
library(CRFutil)


# Cube field model:
grphf <-
  ~X.1:X.2   + X.1:X.4   + X.1:X.5  + X.1:X.5 +
   X.2:X.3 + X.2:X.6 + 
   X.3:X.4 + X.3:X.7 +
   X.4:X.1 + X.4:X.8 +
   X.5:X.6 + X.5:X.8 +
   X.6:X.7 +
   X.7:X.8

adj <- ug(grphf, result="matrix")
adj
node2nme <- data.frame(1:ncol(adj),colnames(adj))
colnames(node2nme) <- c("node","name")
node2nme

# Make up random potentials/sample and return a CRF-object
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}
cub <- make.empty.field(adj.mat = adj, parameterization.typ = "standard", plotQ = T)
#set.seed(87)
cub$par <- runif(cub$n.par,-1.5,1.5)
cub$par # "true" theta
out.pot <- make.pots(parms = cub$par,  crf = cub,  rescaleQ = F, replaceQ = T)
cub$edges
cub$node.pot
cub$edge.pot

# So now sample from the model as if we obtained an experimental sample:
num.samps <- 500
#set.seed(1)
samps <- sample.exact(cub, num.samps)
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

# Model Matrix with respect to graph ????
M <- compute.model.matrix(configs=X.all, edges.mat=cub$edges, node.par = cub$node.par, edge.par = cub$edge.par, ff = f0)
#M
dim(M)

# dat <- list(
#   N = nrow(X.all), # num obs
#   K = ncol(M),     # num vars
#   x = M,           # vars
#   y = freqs        # counts
# )
dat <- list(
  p = ncol(M),
  N = nrow(X.all),
  y = freqs,      
  Mmodl = M       
)
dat

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
#fpth <- "C:/Users/aliso/codes/CRFutil/tests/regression_tests/hojsgaard_model_tests/zero_inflated_fish/"
#fpth <- "/home/npetraco/codes/R/CRFutil/tests/regression_tests/hojsgaard_model_tests/zero_inflated_fish/"
#model.c <- stanc(file = paste0(fpth,"zero_inflated_poisson.stan"), model_name = 'model')
#model.c <- stanc(file = paste0(fpth,"zero_inflated_poisson.stan"), model_name = 'model')
#fpth <- "C:/Users/aliso/codes/CRFutil/tests/regression_tests/hojsgaard_model_tests/triangle_model/"
fpth <- "/home/npetraco/codes/R/CRFutil/tests/regression_tests/hojsgaard_model_tests/triangle_model/"
model.c <- stanc(file = paste0(fpth,"vanalla.poisson.regression2a.stan"), model_name = 'model')
sm <- stan_model(stanc_ret = model.c, verbose = T)

#fit <- sampling(sm, data = dat)
#, control=list(adapt_delta=0.95, max_treedepth=20)
fit <- sampling(sm, data = dat, chains = 4, iter = 2000, thin=1)
#options(max.print = 9999999)
fit.smy <- summary(fit)$summary
rhats <- fit.smy[,"Rhat"]
neffs <- fit.smy[,"n_eff"]

