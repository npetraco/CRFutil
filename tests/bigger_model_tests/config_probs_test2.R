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
fit <- fit_mle_params(grphf, samps,
                      parameterization.typ = "standard",
                      opt.method           = "L-BFGS-B",
                      inference.method     = infer.trbp,
                      state.nmes           = c("white","black"),
                      num.iter             = 5,
                      mag.grad.tol         = 1e-3)
class(fit)
fit$node.potentials
fit$edge.potentials
fit$node.energies
fit$edge.energies
grphf
fit


logZ <- infer.trbp(fit)$logZ
logZ
# Make up a configuration:
X  <- sample(x = c("bk","wt"), size = 6*8, replace = T)

# \Pr({\bf X}) = \frac{1}{Z} e^{E({\bf X})}
EX <- config.energy(config       = X,
                    edges.mat    = fit$edges,
                    one.lgp      = fit$node.energies,
                    two.lgp      = fit$edge.energies, # make sure use same order as edges!
                    ff           = f0)
EX - logZ       # log(Pr(X))
exp(EX - logZ)  # Pr(X)


# Most likely config?
data(Small)
d <- decode.trbp(Small$crf, verbose = T)
d
