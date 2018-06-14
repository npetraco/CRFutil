library(CRFutil)

# Slayer field:
grphf <- ~A:B+A:C+A:D+A:E+B:C+B:D+B:E+C:D+D:E
gp    <- ug(grphf, result = "graph")
adj   <- ug(grphf, result="matrix")

#Make up random potentials and return a CRF-object and sample from it and fit to get a theta:
slay  <- sim.field.random(adjacentcy.matrix=adj, num.states=2, num.sims=100, seed=1)$model
samps <- sample.exact(slay, size = 1000)
fit   <- mrf.standard.fit(samps, grphf, num.states=2, mrf.exact.nll, infer.exact)$fit.model
theta <- fit$par
theta

# Define states and feature function:
s1 <- 1
s2 <- 2
f  <- function(y){ as.numeric(c((y==s1),(y==s2))) }

# Compute features for configuration:
config.phi.vec <- phi.features(config = samps[87, ], fit$edges, fit$node.par, fit$edge.par, f)

theta %*% config.phi.vec

# Decorate parameters
log(fit$node.pot)
fit.pe$node.energies

fit.pe <- make.gRbase.potentials(crf=fit, node.names=gp@nodes, state.nmes=c(s1,s2))
config.energy(config = samps[87, ], fit$edges, fit.pe$node.energies, fit.pe$edge.energies, f)

theta2 <- theta
theta2[c(3,4)] <- -1*theta2[c(3,4)]
theta2
theta2 %*% config.phi.vec

unlist(fit.pe$node.energies)
theta
fit.pe$edge.energies

M <- compute.model.matrix(samps, fit$edges, fit$node.par, fit$edge.par, f)
colSums(M)
fit$par.stat
