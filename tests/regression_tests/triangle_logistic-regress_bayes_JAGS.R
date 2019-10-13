library(CRFutil)
library(CRFutilRcppComponents)
library(R2jags)

# Triangle model
grphf <- ~1:2+1:3+2:3
adj <- ug(grphf, result="matrix")
gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)

f0 <- function(y){ as.numeric(c((y==1),(y==2)))}
n.states <- 2

# Instantiate model and pick a "true" theta:
knm <- make.empty.field(adj.mat = adj, parameterization.typ = "standard", plotQ = T)
set.seed(6)
knm$par <- runif(6,-1.5,1.1)
knm$par
# "true" theta:
#0.07629757  0.93786913 -0.81268462 -0.51175581  0.59945681  1.04299689
out.pot <- make.pots(parms = knm$par,  crf = knm,  rescaleQ = T, replaceQ = T)
# So now sample from the model as if we obtained an experimental sample:
num.samps <- 25
set.seed(1)
samps <- sample.exact(knm, num.samps)
mrf.sample.plot(samps)

# Prepare a logistic model to fit. This time we will use JAGS instead of Stan
fit <- make.empty.field(adj.mat = adj, parameterization.typ = "standard", plotQ = T)

# Here we remove the extra index put in by CRF using the C function. We need to do this for the the RcppComponent functions
theta.pars   <- fix_node_and_edge_par(node_par = fit$node.par, edge_par = fit$edge.par)
fit$edge.par # Should we replace with what's in theta.pars??

Da.mat <- delta_alpha(
  samples       = as.matrix(samps),
  node_par      = theta.pars$node_par,
  edge_par      = theta.pars$edge_par,
  edge_mat      = fit$edges,
  adj_nodes     = fit$adj.nodes,
  num_params_in = 0)

Da.mat

# Stack the node responses into one long binary vector:
y <-c(samps[,1], samps[,2], samps[,3])
y[which(y==2)] <- 0
y

K <- ncol(Da.mat)
N <- nrow(Da.mat)
count.data <- list (y=y, Delta_alpha=Da.mat, N=N, K=K)
count.data
names(count.data)
count.parameters <- c ("beta")

count.inits <- function (){
  #list(beta=rnorm(degree), offset=rnorm(1), sig.beta=runif(1), eps=rnorm(K))
  list(theta=rnorm(K))
  #  list(beta=rnorm(degree), offset=rnorm(1), sig.beta=runif(1), eps=rnorm(K)),
  #  list(beta=rnorm(degree), offset=rnorm(1), sig.beta=runif(1), eps=rnorm(K))
}

model.file.name <- file.path("inst/logistic_model.bug")

#Run model with JAGS instead:
jsim <- jags(data = count.data,
             inits = count.inits,
             parameters.to.save = count.parameters,
             n.iter=10, n.chains=1, model.file=model.file.name)
