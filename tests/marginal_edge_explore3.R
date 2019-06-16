library(CRFutil)

# Using the same node potentials (0.8,0.2), (0.5,0.5), what does the re-scaled edge potential look like under independence??
num.samps <- 100000

# Must label sample states 1 and 2
samps <- cbind(
  1+rbinom(n = num.samps, size = 1, prob = 0.2),
  1+rbinom(n = num.samps, size = 1, prob = 0.5)
)


mrf.sample.plot(samps)
samps

junk <- marginal.edge.info(samps)
dump.crf(crf = junk)
class(samps)

# MLE by hand
junk.mrf <- make.empty.field(graph.eq = ~X1:X2, parameterization.typ = "standard", plotQ = F)
#junk.mrf <-train.mrf(junk.mrf, samps)
#samps

# Auxiliary, gradient convenience function. Follows train.mrf in CRF:
gradient <- function(par, crf, ...) { crf$gradient }
junk.mrf$par.stat <- mrf.stat(junk.mrf, instances = samps)
#junk.mrf$par.stat

# Lets use loopy-belief (lbp) to compute any needed inference quantities (Z and Bels)
# I had to run optim 3-times to reach convergence with LBP:
infr.meth <- infer.exact        # inference method needed for Z and marginals calcs
opt.info  <- stats::optim(    # optimize parameters
  par          = junk.mrf$par,       # theta
  fn           = negloglik,     # objective function
  gr           = gradient,      # grad of obj func
  crf          = junk.mrf,           # passed to fn/gr
  samples      = samps,         # passed to fn/gr
  infer.method = infr.meth,     # passed to fn/gr
  update.crfQ  = TRUE,          # passed to fn/gr
  method       = "L-BFGS-B",
  control      = list(factr=10, pgtol=1e-16, trace = 1, REPORT=1))
opt.info$convergence
opt.info$message
junk.mrf$gradient
junk.mrf$nll

dump.crf(junk.mrf)
