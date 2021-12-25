library(CRFutil)
library(plot.matrix)

# Graph: THIS IS TOO BIG TO COMPUTE THE FULL JOINT........ 2^36 states MUST BE ENUMERATED.... Can get Z though
grphf      <- make.lattice(num.rows = 6, num.cols = 6, cross.linksQ = T)
adj        <- ug(grphf, result="matrix") # adjacency (connection) matrix
node.names <- colnames(adj)
# Check the graph:
gp <- ug(grphf, result = "graph")
plot(gp)
dev.off()

# Make up random parameters for the graph and simulate some data from it:
known.model.info <- sim.field.random(adjacentcy.matrix=adj, num.states=2, num.sims=25, seed=1)
samps            <- known.model.info$samples
known.model      <- known.model.info$model

# Fit an MRF to the sample with the intention of obtaining the estimated parameter vector theta
# Use the standard parameterization (one parameter per node, one parameter per edge):
s1 <- "white" # State 1 name
s2 <- "black" # State 2 name
# s1 <- 1 # State 1 name
# s2 <- 2 # State 2 name
fit <- fit_mle_params(grphf, samps,
                      parameterization.typ = "standard",
                      opt.method           = "L-BFGS-B",
                      #opt.method           = "CG",
                      inference.method     = infer.trbp,
                      state.nmes           = c(s1,s2),
                      num.iter             = 1,
                      mag.grad.tol         = 1e-3,
                      plotQ                = T)

# Prep for computing config probs:
logZ <- infer.trbp(fit)$logZ
logZ
f0   <- function(y){ as.numeric(c((y==s1),(y==s2))) }

# Make up a configuration:
X <- sample(x = c(s1,s2), size = fit$n.nodes, replace = T)
# See what the configuration looks like
Xmat <- t(array(X, c(sqrt(fit$n.nodes),sqrt(fit$n.nodes))))
Xmat
par(mar=c(5.1, 4.1, 4.1, 4.1))
plot(Xmat, col=c("white", "black"))


# Compute the configuration probability: \Pr({\bf X}) = \frac{1}{Z} e^{E({\bf X})}
EX <- config.energy(config       = X,
                    edges.mat    = fit$edges,
                    one.lgp      = fit$node.energies,
                    two.lgp      = fit$edge.energies, # make sure use same order as edges!
                    ff           = f0)
EX
EX - logZ       # log(Pr(X))
exp(EX - logZ)  # Pr(X)


# Try some of the sampled configurations. They probably have higher Pr(X)'s:
X <- samps[1,]
X[which(X == 1)] <- s1
X[which(X == 2)] <- s2
X
# See what the configuration looks like
Xmat <- t(array(X, c(sqrt(fit$n.nodes),sqrt(fit$n.nodes))))
Xmat
plot(Xmat, col=c("white", "black"))

# Pr(Xmat):
EX <- config.energy(config       = X,
                    edges.mat    = fit$edges,
                    one.lgp      = fit$node.energies,
                    two.lgp      = fit$edge.energies, # make sure use same order as edges!
                    ff           = f0)
EX
EX - logZ       # log(Pr(X))
exp(EX - logZ)  # Pr(X)


# Most likely config?
#decode.trbp(fit, verbose = T)
decode.junction(fit)
decode.greedy(fit)
#decode.exact(fit)
#decode.tree(fit)
decode.lbp(fit)


# Pr of most likely configuration??:
X <- decode.junction(fit)
X
X[which(X == 1)] <- s1
X[which(X == 2)] <- s2
X
# See what the configuration looks like
Xmat <- t(array(X, c(sqrt(fit$n.nodes),sqrt(fit$n.nodes))))
Xmat
plot(Xmat, col=c("white", "black"))

# Pr(Most likely config): Should be relatively high compared to other configs
EX <- config.energy(config       = X,
                    edges.mat    = fit$edges,
                    one.lgp      = fit$node.energies,
                    two.lgp      = fit$edge.energies, # make sure use same order as edges!
                    ff           = f0)
EX
EX - logZ       # log(Pr(X))
exp(EX - logZ)  # Pr(X)


# This model is (just) small enough that we can examine the estimated joint distribution and compare
# it with the true joint distribution:
# Compute the joint dist from the estimated theta:
gR.dist.info.mle.model <- distribution.calc(fit, logZ.calc = infer.junction)
gR.dist.info.mle.model$logZ

logZ.mle.model         <- gR.dist.info.mle.model$logZ # Uses exact computation for logZ
logZ.mle.model

# Just to see if this works:
infer.exact(fit)$logZ

joint.mle <- as.data.frame(as.table(gR.dist.info.mle.model$state.probs))
dim(joint.mle)

pr.idx <- ncol(joint.mle)
#round(joint.mle[,pr.idx],3)
#plot(joint.mle[,pr.idx], typ="h")

# TRUE joint distribution:
# Compare estimated joint distribution to the true distribution
# Now compute the joint dist from the estimated theta
out.pot.true <- make.pots(parms = known.model$par,  crf = known.model,  rescaleQ = F, replaceQ = T)
gR.true      <- make.gRbase.potentials(known.model, node.names = node.names, state.nmes = c(s1,s2))

known.model$node.potentials <- gR.true$node.potentials
known.model$edge.potentials <- gR.true$edge.potentials

gR.dist.info.true.model <- distribution.calc(known.model, logZ.calc = infer.exact)

logZ.true.model <- gR.dist.info.true.model$logZ
logZ.true.model

joint.true      <- as.data.frame(as.table(gR.dist.info.true.model$state.probs))

colnames(joint.true)
colnames(joint.mle)    # NODE ORDERING THE SAME ??????

# cbind(
#   round(joint.mle[,pr.idx],3),
#   round(joint.true[,pr.idx],3)
# )

# Distribution discrepancy measures
# par(mfrow=c(2,1))
# ymax <- max(joint.true[,pr.idx], joint.mle[,pr.idx])
# plot(joint.true[,pr.idx], typ="h", ylim=c(0,ymax))
# plot(joint.mle[,pr.idx],  typ="h", ylim=c(0,ymax))
# dev.off()

# plot(abs(joint.true[,pr.idx] - joint.mle[,pr.idx]), typ="h", ylab="abs-diff (Pr-units)")
# plot(joint.true[,pr.idx] - joint.mle[,pr.idx], typ="h", ylab="diff (Pr-units)")
#plot(abs(joint.true[,pr.idx] - joint.mle[,pr.idx])/joint.true[,pr.idx]*100, typ="h", ylab="%-diff (%)")

# Kullback-Leibler distance between true and estimated distribution
KLD(joint.true[,pr.idx], joint.mle[,pr.idx])

# Jensen-Shannon distance between true and estimated distribution (Should be bound between 0 and 1)
PD  <- joint.true[,pr.idx]
QD  <- joint.mle[,pr.idx]
RD  <- 0.5*(joint.true[,pr.idx] + joint.mle[,pr.idx])
JSD <- 0.5 * (KLD(PD,RD,base = 2)$sum.KLD.px.py + KLD(QD,RD, base = 2)$sum.KLD.px.py)
JSD

