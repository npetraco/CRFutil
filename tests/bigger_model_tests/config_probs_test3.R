library(CRFutil)

# Graph:
grphf      <- make.lattice(num.rows = 6, num.cols = 8, cross.linksQ = T)
adj        <- ug(grphf, result="matrix") # adjacency (connection) matrix
node.names <- colnames(adj)
# Check the graph:
gp <- ug(grphf, result = "graph")
plot(gp)
dev.off()

# Make up random parameters for the graph and simulate some data from it:
known.model.info <- sim.field.random(adjacentcy.matrix=adj, num.states=2, num.sims=25)
samps            <- known.model.info$samples
known.model      <- known.model.info$model

# Fit an MRF to the sample with the intention of obtaining the estimated parameter vector theta
# Use the standard parameterization (one parameter per node, one parameter per edge):
s1 <- "white" # State 1 name
s2 <- "black" # State 2 name
fit <- fit_mle_params(grphf, samps,
                      parameterization.typ = "standard",
                      opt.method           = "L-BFGS-B",
                      #opt.method           = "CG",
                      inference.method     = infer.trbp,
                      state.nmes           = c(s1,s2),
                      num.iter             = 10,
                      mag.grad.tol         = 1e-3,
                      plotQ                = T)

# Prep for computing config probs:
logZ <- infer.trbp(fit)$logZ
logZ
f0   <- function(y){ as.numeric(c((y==s1),(y==s2))) }

# Make up a configuration:
X <- sample(x = c(s1,s2), size = 6*8, replace = T)

# \Pr({\bf X}) = \frac{1}{Z} e^{E({\bf X})}
EX <- config.energy(config       = X,
                    edges.mat    = fit$edges,
                    one.lgp      = fit$node.energies,
                    two.lgp      = fit$edge.energies, # make sure use same order as edges!
                    ff           = f0)
EX
EX - logZ       # log(Pr(X))
exp(EX - logZ)  # Pr(X)


# Try some of the sampled configurations:
X <- samps[1,]
X[which(X == 1)] <- s1
X[which(X == 2)] <- s2
X

EX <- config.energy(config       = X,
                    edges.mat    = fit$edges,
                    one.lgp      = fit$node.energies,
                    two.lgp      = fit$edge.energies, # make sure use same order as edges!
                    ff           = f0)
EX
EX - logZ       # log(Pr(X))
exp(EX - logZ)  # Pr(X)


# Most likely config?
decode.trbp(fit, verbose = T)
decode.junction(fit)
decode.greedy(fit)
decode.tree(fit)
decode.lbp(fit)

# Pr of most likely configuration??:
X <- decode.junction(fit)
X
X[which(X == 1)] <- s1
X[which(X == 2)] <- s2
X
EX <- config.energy(config       = X,
                    edges.mat    = fit$edges,
                    one.lgp      = fit$node.energies,
                    two.lgp      = fit$edge.energies, # make sure use same order as edges!
                    ff           = f0)
EX
EX - logZ       # log(Pr(X))
exp(EX - logZ)  # Pr(X)
