library(CRFutil)

knm <- sim.random.field(
  num.nodes=10, num.edges=13, plotQ = T,
  node.pot.nu   = 30,
  node.pot.mean = 1,
  node.pot.sd   = 3,
  edge.pot.nu   = 30,
  edge.pot.mean = 0,
  edge.pot.sd   = 3,
  num.sims      = 100)
knm.crf <- knm$model
samps   <- knm$samples
mrf.sample.plot(samps)
plot_crf(knm.crf)

knm.dist <- compute.full.distribution(knm.crf)$joint.distribution


# Start with saturated graph
g.start             <- erdos.renyi.game(ncol(samps), 1, typ="gnp")
adj.start           <- as.matrix(as_adj(g.start))
colnames(adj.start) <- 1:ncol(samps)
rownames(adj.start) <- 1:ncol(samps)
satm                 <- make.empty.field(adj.mat = adj.start, parameterization.typ = "standard", plotQ = T)

# Compute the delta-alpha matrix given the sampled configs:   # ***********ADD nfolds arguement
f0               <- function(y){ as.numeric(c((y==1),(y==2))) }
fit.glmn.logis   <- mrf.glmnet.logistic.fit(samples = samps, crf.obj = satm, lambda="lambda.1se", ff=f0, infer.func=infer.junction, plotQ=T)
#

# MLE fit for parameters of glmnet graph:
fit.glmn.logis.mle.adj <- edges2adj(fit.glmn.logis$edges, n.nodes = fit.glmn.logis$n.nodes, plotQ = F)
fit.glmn.logis.mle <- make.empty.field(adj.mat = fit.glmn.logis.mle.adj, parameterization.typ = "standard", plotQ = T)

# Compute the sufficient stats needed by the likelihood and its’ grad
fit.glmn.logis.mle$par.stat <- mrf.stat(fit.glmn.logis.mle, samps)

# Auxiliary, gradient convenience function. Follows train.mrf in CRF:
gradient <- function(par, crf, ...) { crf$gradient }

# MLE on the specified graph given the samples:
infr.meth <- infer.exact   # inference method needed for Z and marginals calcs
opt.info  <- stats::optim(    # optimize parameters
  par          = fit.glmn.logis.mle$par,       # theta
  fn           = negloglik,     # objective function
  gr           = gradient,      # grad of obj func
  crf          = fit.glmn.logis.mle,           # passed to fn/gr
  samples      = samps,         # passed to fn/gr
  infer.method = infr.meth,     # passed to fn/gr
  update.crfQ  = TRUE,          # passed to fn/gr
  method       = "L-BFGS-B",
  control      = list(trace = 1, REPORT=1))
opt.info$convergence
opt.info$message
fit.glmn.logis.mle$gradient
fit.glmn.logis.mle$nll

dump.crf(fit.glmn.logis.mle)

fit.glmn.logis.dist <- compute.full.distribution(fit.glmn.logis)$joint.distribution
fit.glmn.logis.mle.dist <- compute.full.distribution(fit.glmn.logis.mle)$joint.distribution

plot(1:nrow(knm.dist), knm.dist$Freq, typ="h", col="green")           # True joint dist
points((1:nrow(fit.glmn.logis.dist))+0.5, fit.glmn.logis.dist$Freq, typ="h", col="red") # glmnet selected joint dist
points((1:nrow(fit.glmn.logis.mle.dist))-0.5, fit.glmn.logis.mle.dist$Freq, typ="h", col="blue") # glmnet selected joint dist

kdist.info <- KLD(knm.dist$Freq, fit.glmn.logis.dist$Freq)
kdist.info$sum.KLD.px.py
kdist.info$sum.KLD.py.px
kdist.info$mean.sum.KLD
kdist.info$intrinsic.discrepancy

kdist.info2 <- KLD(knm.dist$Freq, fit.glmn.logis.mle.dist$Freq)
kdist.info2$sum.KLD.px.py
kdist.info2$sum.KLD.py.px
kdist.info2$mean.sum.KLD
kdist.info2$intrinsic.discrepancy

compare_edges(model.list=list(knm.crf, fit.glmn.logis), num.nodes=NULL)
plot_crf(knm.crf)
plot_crf(fit.glmn.logis)
fit.glmn.logis$edges
