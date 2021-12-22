library(CRFutil)

# Graph:
grphf      <- make.lattice(num.rows = 32, num.cols = 32, cross.linksQ = T)
adj        <- ug(grphf, result="matrix") # adjacency (connection) matrix
node.names <- colnames(adj)
# Check the graph:
gp <- ug(grphf, result = "graph")
dev.off()
plot(gp)
dev.off()

# Make up random parameters for the graph and simulate some data from it:
known.model.info <- sim.field.random(adjacentcy.matrix=adj, num.states=2, num.sims=25) #seed=1 for reproducible results
samps            <- known.model.info$samples
known.model      <- known.model.info$model

# Fit an MRF to the sample with the intention of obtaining the estimated parameter vector theta
# Use the standard parameterization (one parameter per node, one parameter per edge):
fit <- make.empty.field(adj.mat = adj, parameterization.typ = "standard", plotQ = F)

# MLE for parameters of model. Follows train.mrf in CRF, just gives more control and output:
gradient      <- function(par, crf, ...) { crf$gradient } # Auxiliary, gradient convenience function.
fit$par.stat <- mrf.stat(fit, samps)   # requisite sufficient statistics
#infr.meth <- infer.exact               # inference method needed for Z and marginals calcs
#infr.meth <- infer.junction            # inference method needed for Z and marginals calcs
#infr.meth <- infer.cutset            # inference method needed for Z and marginals calcs
infr.meth <- infer.trbp            # inference method needed for Z and marginals calcs
#infr.meth <- infer.lbp
#infr.meth <- infer.rbp
#infr.meth <- infer.sample
#infr.meth <- infer.tree
opt.info  <- stats::optim(             # optimize parameters
  par          = fit$par,              # theta
  fn           = negloglik,            # objective function
  gr           = gradient,             # grad of obj func
  crf          = fit,                  # passed to fn/gr
  samples      = samps,                # passed to fn/gr
  infer.method = infr.meth,            # passed to fn/gr
  update.crfQ  = TRUE,                 # passed to fn/gr
  method       = "L-BFGS-B",
  control      = list(pgtol = 1e-4, trace = 1, REPORT=1))

# Checks: May have to re-run optimize a few times to get gradient down:
opt.info$convergence
opt.info$message
fit$gradient
fit$nll
fit$par         # Estimated parameter vector
known.model$par # True parameter vector

# Now compute the joint dist from the estimated theta
out.pot <- make.pots(parms = fit$par,  crf = fit,  rescaleQ = F, replaceQ = T)
gR.mle  <- make.gRbase.potentials(fit, node.names = node.names, state.nmes = c("bk","wt"))

gR.dist.info.mle.model <- distribution.from.potentials(gR.mle$node.potentials,
                                                       gR.mle$edge.potentials)
logZ.mle.model         <- gR.dist.info.mle.model$logZ
joint.mle <- as.data.frame(as.table(gR.dist.info.mle.model$state.probs))

dim(joint.mle)
pr.idx <- ncol(joint.mle)
#round(joint.mle[,pr.idx],3)
plot(joint.mle[,pr.idx], typ="h")

# Compare estimated joint distribution to the true distribution
# Now compute the joint dist from the estimated theta
out.pot.true <- make.pots(parms = known.model$par,  crf = known.model,  rescaleQ = F, replaceQ = T)
gR.true      <- make.gRbase.potentials(known.model, node.names = node.names, state.nmes = c("bk","wt"))

gR.dist.info.true.model <- distribution.from.potentials(gR.true$node.potentials,
                                                        gR.true$edge.potentials)
logZ.true.model <- gR.dist.info.true.model$logZ
joint.true      <- as.data.frame(as.table(gR.dist.info.true.model$state.probs))

#cbind(joint.true, joint.mle)

colnames(joint.true)
colnames(joint.mle)    # NODE ORDERING THE SAME ??????

# cbind(
#   round(joint.mle[,pr.idx],3),
#   round(joint.true[,pr.idx],3)
# )

# Distribution discrepancy measures
par(mfrow=c(2,1))
ymax <- max(joint.true[,pr.idx], joint.mle[,pr.idx])
plot(joint.true[,pr.idx], typ="h", ylim=c(0,ymax))
plot(joint.mle[,pr.idx],  typ="h", ylim=c(0,ymax))
dev.off()

plot(abs(joint.true[,pr.idx] - joint.mle[,pr.idx]), typ="h", ylab="abs-diff (Pr-units)")
plot(joint.true[,pr.idx] - joint.mle[,pr.idx], typ="h", ylab="diff (Pr-units)")
#plot(abs(joint.true[,pr.idx] - joint.mle[,pr.idx])/joint.true[,pr.idx]*100, typ="h", ylab="%-diff (%)")

# Kullback-Leibler distance between true and estimated distribution
KLD(joint.true[,pr.idx], joint.mle[,pr.idx])

# Jensen-Shannon distance between true and estimated distribution (Should be bound between 0 and 1)
PD  <- joint.true[,pr.idx]
QD  <- joint.mle[,pr.idx]
RD  <- 0.5*(joint.true[,pr.idx] + joint.mle[,pr.idx])
JSD <- 0.5 * (KLD(PD,RD,base = 2)$sum.KLD.px.py + KLD(QD,RD, base = 2)$sum.KLD.px.py)
JSD

