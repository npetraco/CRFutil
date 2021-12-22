library(CRFutil)

# Graph:
grphf      <- make.lattice(num.rows = 6, num.cols = 8, cross.linksQ = T)
adj        <- ug(grphf, result="matrix") # adjacency (connection) matrix
node.names <- colnames(adj)
# Check the graph:
gp <- ug(grphf, result = "graph")
dev.off()
plot(gp)
dev.off()

# Make up random parameters for the graph and simulate some data from it:
known.model.info <- sim.field.random(adjacentcy.matrix=adj, num.states=2, num.sims=25)
samps            <- known.model.info$samples
known.model      <- known.model.info$model

# Fit an MRF to the sample with the intention of obtaining the estimated parameter vector theta
# Use the standard parameterization (one parameter per node, one parameter per edge):
fit <- make.empty.field(adj.mat = adj, parameterization.typ = "standard", plotQ = F)

a.list <- list("a","b")
a.list
a.list[[3]] <- "c"
a.list

b.list <- list("d","e")
c(a.list, b.list)

#MLE for parameters of model. Follows train.mrf in CRF, just gives more control and output:
gradient      <- function(par, crf, ...) { crf$gradient } # Auxiliary, gradient convenience function.
fit$par.stat <- mrf.stat(fit, samps)   # requisite sufficient statistics
#infr.meth <- infer.exact              # inference method needed for Z and marginals calcs
#infr.meth <- infer.junction           # inference method needed for Z and marginals calcs
infr.meth <- infer.trbp                # inference method needed for Z and marginals calcs
num.iters <- 15
for(i in 1:num.iters){
  opt.info  <- stats::optim(             # optimize parameters
    par          = fit$par,              # theta
    fn           = negloglik,            # objective function
    gr           = gradient,             # grad of obj func
    crf          = fit,                  # passed to fn/gr
    samples      = samps,                # passed to fn/gr
    infer.method = infr.meth,            # passed to fn/gr
    update.crfQ  = TRUE,                 # passed to fn/gr
    method       = "L-BFGS-B",
    control      = list(trace = 1, REPORT=1))
}

# Checks: May have to re-run optimize a few times to get gradient down:
opt.info$convergence
opt.info$message
fit$gradient
fit$nll
fit$par         # Estimated parameter vector
known.model$par # True parameter vector

# Prep to compute configurtion energies
out.pot <- make.pots(parms = fit$par,  crf = fit,  rescaleQ = F, replaceQ = T)
gR.mle  <- make.gRbase.potentials(fit, node.names = node.names, state.nmes = c("bk","wt"))
#names(gR.mle)
logZ    <- infer.trbp(fit)$logZ

f0      <- function(y){ as.numeric(c((y=="bk"),(y=="wt"))) }

# A configuration:
X  <- sample(x = c("bk","wt"), size = 6*8, replace = T)

# \Pr({\bf X}) = \frac{1}{Z} e^{E({\bf X})}
EX <- config.energy(config = X,
              edges.mat    = fit$edges,
              one.lgp      = gR.mle$node.energies,
              two.lgp      = gR.mle$edge.energies, # make sure use same order as edges!
              ff           = f0)
EX - logZ       # log(Pr(X))
exp(EX - logZ)  # Pr(X)

# Most likely config?
# Plot configs as array
