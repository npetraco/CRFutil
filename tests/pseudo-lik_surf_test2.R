library(CRFutil)
library(rgl)

#
grphf <- ~A:B
adj   <- ug(grphf, result="matrix")
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}

n.states <- 2
known.model <- make.crf(adj, n.states)

# True node pots:
PsiA <- c(1,1)
PsiB <- c(1,1)

# True edge pots:
PsiAB <-
  rbind(
    c(5, 2.5),
    c(2.5, 5)
  )

known.model$node.pot[1,]  <- PsiA
known.model$node.pot[2,]  <- PsiB
known.model$edge.pot[[1]] <- PsiAB

# So now sample from the model as if we obtained an experimental sample:
num.samps <- 25
set.seed(1)
samps <- sample.exact(known.model, num.samps)
samps
mrf.sample.plot(samps)

fit <- make.crf(adj, n.states)
fit <- make.features(fit)
fit <- make.par(fit, 2)

fit$node.par[1,1,] <- 1
fit$node.par[2,1,] <- 1
fit$node.par

# Use the more flexible edge parameterization:
fit$edge.par[[1]][1,1,1] <- 2
fit$edge.par[[1]][2,2,1] <- 2
fit$edge.par

# First compute the sufficient stats needed by the likelihood and its’ grad
fit$par.stat <- mrf.stat(fit, samps)
fit$par.stat

# Auxiliary, gradient convenience function. Follows train.mrf in CRF:
gradient <- function(par, crf, ...) { crf$gradient }

# Lets use loopy-belief (lbp) to compute any needed inference quantities (Z and Bels)
# I had to run optim 3-times to reach convergence with LBP:
infr.meth <- infer.exact        # inference method needed for Z and marginals calcs
opt.info  <- stats::optim(    # optimize parameters
  par          = fit$par,       # theta
  fn           = negloglik,     # objective function
  gr           = gradient,      # grad obj func
  crf          = fit,           # passed to fn/gr
  samples      = samps,         # passed to fn/gr
  infer.method = infr.meth,     # passed to fn/gr
  update.crfQ  = TRUE,          # passed to fn/gr
  method       = "L-BFGS-B",
  control      = list(trace = 1, REPORT=1))
opt.info$convergence
opt.info$message
fit$gradient
fit$nll

infer.exact(fit)
fit$node.pot
fit$edge.pot

fit$par
#----------------------------------------------------

fit2 <- make.crf(adj, n.states)
fit2 <- make.features(fit2)
fit2 <- make.par(fit2, 2)
fit2$node.par[1,1,] <- 1
fit2$node.par[2,1,] <- 1
fit2$node.par

fit2$edge.par[[1]][1,1,1] <- 2
fit2$edge.par[[1]][2,2,1] <- 2
fit2$edge.par

fit$par
theta1 <- seq(from = -1.75, 0.0, by = 0.05)
theta2 <- seq(from = -1.75, 0.0, by = 0.05)
length(theta1)
length(theta1)^2
theta.grid <- as.matrix(expand.grid(theta1,theta2))
dim(theta.grid)

n2p <- nodes2params.list(fit2, storeQ = T)

# First lets try to compute the negative log likelihood for a FIXED parameter vector theta:
theta.fixed <- theta.grid[3055,]
theta.fixed

fit2$par
neg.log.psl <- sum(sapply(1:nrow(samps), function(xx){
  neglogpseudolik.config(par = theta.fixed, config = samps[xx,], crf = fit2, ff = f0)}))
neg.log.psl

# Now lets try looping over the different parameter vectors contained in theta.grid. For each compute the
# corresponding sample pseudolikelihood:
neg.log.psls <- array(NA, nrow(theta.grid))
for(i in 1:nrow(theta.grid)) {

  print(i)
  theta <- theta.grid[i,]
  neg.log.psls[i] <- sum(sapply(1:nrow(samps), function(xx){neglogpseudolik.config(par = theta, config = samps[xx,], crf = fit2, ff = f0)}))

}
neg.log.psls
hist(neg.log.psls)
min(neg.log.psls)
zw <- theta.grid[which(neg.log.psls < 33),]
c(min(zw[,1]), max(zw[,1]))
c(min(zw[,2]), max(zw[,2]))

plot3d(theta.grid[,1], theta.grid[,2], neg.log.psls)

