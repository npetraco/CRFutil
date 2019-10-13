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

#model_data = list(TT = N, y = y, x = Da.mat, K = 1, pp=ncol(Da.mat)) # v1a
model_data = list(TT = N, y = y, x = Da.mat, pp=ncol(Da.mat)) # v1b
model_data

fpth <- "inst/logistic_model_v1b.bug"
model_parameters =  c("beta")

# Run the model
model_run = jags(data = model_data, #inits = init,
                 parameters.to.save = model_parameters,
                 model.file = fpth,
                 n.chains = 4,
                 n.iter = 10000,
                 n.burnin = 1000,
                 n.thin = 10)


#

# Check the output - are the true values inside the 95% CI?
# Also look at the R-hat values - they need to be close to 1 if convergence has been achieved
plot(model_run)
print(model_run)
traceplot(model_run)

# Create a plot of the posterior mean regression line
post = print(model_run)
#beta_1_mean = post$mean$beta_1
#beta_2_mean = post$mean$beta_2
fit$par <- post$mean$beta

out.pot2 <- make.pots(parms = fit$par,  crf = fit,  rescaleQ = T, replaceQ = T)
fit$node.pot
fit$edge.pot

# Configuration probabilities:
pot.info <- make.gRbase.potentials(fit, node.names = gp@nodes, state.nmes = c("1","2"))
pot.info

gR.dist.info    <- distribution.from.potentials(pot.info$node.potentials, pot.info$edge.potentials)
logZ            <- gR.dist.info$logZ
joint.dist.info <- as.data.frame(as.table(gR.dist.info$state.probs))
joint.dist.info

# True model:
pot.info.true.model <- make.gRbase.potentials(knm, node.names = gp@nodes, state.nmes = c("1","2"))
gR.dist.info.true.model    <- distribution.from.potentials(pot.info.true.model$node.potentials, pot.info.true.model$edge.potentials)
logZ.true.model            <- gR.dist.info.true.model$logZ
joint.dist.info.true.model <- as.data.frame(as.table(gR.dist.info.true.model$state.probs))
joint.dist.info.true.model

# Compare:
fit.cp  <- round(joint.dist.info[,4]*100, 1)     # Bayes logistic from JAGS
true.model.cp  <- round(joint.dist.info.true.model[,4]*100, 1) # True
cbind(joint.dist.info[,c(2,3,1)], fit.cp, true.model.cp)

# X1 X2 X3 fit.cp true.model.cp
# 1  1  1  1   10.2          18.4
# 2  1  1  2    1.1           8.0
# 3  2  1  1   22.2          15.6
# 4  2  1  2   23.4          22.6
# 5  1  2  1    8.9           4.2
# 6  1  2  2   10.1          14.9
# 7  2  2  1    2.0           1.3
# 8  2  2  2   22.1          15.0

# Double check
cbind(
  joint.dist.info.true.model,
  joint.dist.info
)
