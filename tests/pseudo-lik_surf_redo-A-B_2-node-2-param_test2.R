library(CRFutil)
library(rgl)

grphf <- ~A:B
adj   <- ug(grphf, result="matrix")
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}

n.states <- 2
known.model <- make.crf(adj, n.states)

# First lets get a sample of configurations from a "true" model:
# The "true" theta is: (1.098612, -0.7096765) which translates to the potentials below.
# True node pots:
PsiA <- c(3,1)
PsiB <- c(3,1)

# True edge pots:
PsiAB <-
  rbind(
    c(3, 6.1),
    c(6.1, 3)
  )

known.model$node.pot[1,]  <- PsiA
known.model$node.pot[2,]  <- PsiB
known.model$edge.pot[[1]] <- PsiAB

# So now sample from the model as if we obtained an experimental sample:
num.samps <- 25
set.seed(1)
samps <- sample.exact(known.model, num.samps)
mrf.sample.plot(samps)

# Now examine the pseudo-likelihood and likelihood surfaces given these samples over a grid of parameter values

# Instantiate a two node two parameter model to fit:
fit <- make.crf(adj, n.states)
fit <- make.features(fit)
fit <- make.par(fit, 2)
fit$node.par[1,1,] <- 1
fit$node.par[2,1,] <- 1
fit$node.par

fit$edge.par[[1]][1,1,1] <- 2
fit$edge.par[[1]][2,2,1] <- 2
fit$edge.par
#n2p <- nodes2params.list2(fit, storeQ = T)


# Set up a grid of theta's overwhich to compute the sample pseudolikelihood:
log(known.model$node.pot)                                            # True theta1
log(known.model$edge.pot[[1]]) - max(log(known.model$edge.pot[[1]])) # True theta2

theta1 <- seq(from =   0.0, to=2.5, by = 0.05)
theta2 <- seq(from =  -2.0, to=1.0, by = 0.05)
theta.grid <- as.matrix(expand.grid(theta1,theta2))
dim(theta.grid)

# Compute the sample log pseudolikelihoods and likelihoods over the theta grid:
neg.log.psls2 <- array(NA, nrow(theta.grid))
neg.log.lik   <- array(NA, nrow(theta.grid))
fit$par.stat <- mrf.stat(crf = fit, instances = samps)

for(i in 1:nrow(theta.grid)) {
  print(i)
  theta.i          <- theta.grid[i,]
  npsls2.info      <- neglogpseudolik(par = theta.i, crf = fit, samples = samps, conditional.energy.function.type = "feature", ff = f0, update.crfQ = F)
  neg.log.psls2[i] <- npsls2.info$samp.neglogpseudolik
  neg.log.lik[i]   <- negloglik(par = theta.i, crf = fit, samples = samps, infer.method = infer.exact, update.crfQ = F)
}

# Look for the rough minimum:
#neg.log.psls2
hist(neg.log.psls2)
min(neg.log.psls2)
max(neg.log.psls2)
min.idx <- which(neg.log.psls2 == min(neg.log.psls2))
# Rough minimum:
c(theta.grid[min.idx,1], theta.grid[min.idx,2], neg.log.psls2[min.idx])
# Minimum in the right ball park?
log(known.model$node.pot)                                            # True theta1
log(known.model$edge.pot[[1]]) - max(log(known.model$edge.pot[[1]])) # True theta2

# Look for the rough minimum of the likelihood formulation:
hist(neg.log.lik)
min(neg.log.lik)
max(neg.log.lik)
min.idx.lik <- which(neg.log.lik == min(neg.log.lik))
# Rough minimum:
c(theta.grid[min.idx.lik,1], theta.grid[min.idx.lik,2], neg.log.lik[min.idx.lik])
# Minimum in the right ball park?
log(known.model$node.pot)                                            # True theta1
log(known.model$edge.pot[[1]]) - max(log(known.model$edge.pot[[1]])) # True theta2


# Plot the sample pseudo-likelihood over the theta grid and the rough minimum:
plot3d(theta.grid[,1], theta.grid[,2], neg.log.psls2, xlab="theta1", ylab="theta2", zlab="", main="feature", main="Neg-log-pseudo-lik")
text3d(theta.grid[min.idx,1], theta.grid[min.idx,2], neg.log.psls2[min.idx], texts = "min")
lines3d(c(theta.grid[min.idx,1],theta.grid[min.idx,1]), c(theta.grid[min.idx,2],theta.grid[min.idx,2]), c(neg.log.psls2[min.idx], max(neg.log.psls2)))
text3d(theta.grid[min.idx,1], theta.grid[min.idx,2], max(neg.log.psls2), texts = "min")

# Examine the likelihood now:
open3d()
# Plot the sample likelihood over the theta grid and the rough minimum:
plot3d(theta.grid[,1], theta.grid[,2], neg.log.lik, xlab="theta1", ylab="theta2", zlab="", main="Neg-log-lik")
text3d(theta.grid[min.idx,1], theta.grid[min.idx,2], neg.log.lik[min.idx], texts = "min")
lines3d(c(theta.grid[min.idx,1],theta.grid[min.idx,1]), c(theta.grid[min.idx,2],theta.grid[min.idx,2]), c(neg.log.lik[min.idx], max(neg.log.lik)))
text3d(theta.grid[min.idx,1], theta.grid[min.idx,2], max(neg.log.lik), texts = "min")


# Gradient of the sample neg log pseudo-likelihood and neg log likelihood:
# Compute and plot gradient:
npll.mat        <- array(NA, c(length(theta1), length(theta2)))
grad.npll.mat.u <- array(NA, c(length(theta1), length(theta2)))
grad.npll.mat.v <- array(NA, c(length(theta1), length(theta2)))
#
nll.mat        <- array(NA, c(length(theta1), length(theta2)))
grad.nll.mat.u <- array(NA, c(length(theta1), length(theta2)))
grad.nll.mat.v <- array(NA, c(length(theta1), length(theta2)))
#
CRF.nll.mat        <- array(NA, c(length(theta1), length(theta2)))
CRF.grad.nll.mat.u <- array(NA, c(length(theta1), length(theta2)))
CRF.grad.nll.mat.v <- array(NA, c(length(theta1), length(theta2)))


count <- 1
for(i in 1:length(theta1)) {
  for(j in 1:length(theta2)) {


    print(count)
    theta            <- c(theta1[i], theta2[j])

    # Negative log pseudo likelihood and gradient:
    npll.info        <- neglogpseudolik(par = theta, crf = fit, samples = samps, conditional.energy.function.type = "feature", ff = f0, gradQ = T, update.crfQ = F)
    npll.mat[i,j]    <- npll.info$samp.neglogpseudolik
    neg.log.psl.grad <- npll.info$samp.grad.neglogpseudolik

    grad.npll.mat.u[i,j] <- neg.log.psl.grad[1]
    grad.npll.mat.v[i,j] <- neg.log.psl.grad[2]

    # Negative log likelihood and gradient:
    # update.crf=T Must be true to use the current theta. NOTE: par in fit will ne changing!
    nll.mat[i,j]        <- negloglik(par = theta, crf = fit, samples = samps, infer.method = infer.exact, update.crfQ = T)
    #neg.log.lik.grad    <- grad.negloglik(crf = fit, nInstances = samps, suffStat = fit$par.stat, inference.func = infer.exact)
    #grad.nll.mat.u[i,j] <- neg.log.lik.grad[1]
    #grad.nll.mat.v[i,j] <- neg.log.lik.grad[2]
    grad.nll.mat.u[i,j] <- fit$gradient[1]
    grad.nll.mat.v[i,j] <- fit$gradient[2]


    # CHECK: nll and gnll from the CRF package:
    #CRF.nll.mat[i,j]        <- mrf.nll(par = theta, crf = fit, instances = samps, infer.method = infer.exact)
    #CRF.grad.nll.mat.u[i,j] <- fit$gradient[1]
    #CRF.grad.nll.mat.v[i,j] <- fit$gradient[2]

    count <- count + 1
  }
}



# Make a 2D plot of the gradient for the negative log pseudo likelihood:
library(OceanView)
u <- grad.npll.mat.u
v <- grad.npll.mat.v

# Sample negative log pseudo-likelihood:
image2D(x = theta1, y = theta2, z = npll.mat, contour = TRUE, xlab="theta1", ylab="theta2")

# Gradient:
quiver2D(x = theta1, y = theta2, u = u, v = v, add = TRUE, by = 6)

#-------
# Make a 2D plot of the gradient for the negative log likelihood:
ul <- grad.nll.mat.u
vl <- grad.nll.mat.v

# Sample negative log likelihood:
image2D(x = theta1, y = theta2, z = nll.mat, contour = TRUE, xlab="theta1", ylab="theta2")

# Gradient:
quiver2D(x = theta1, y = theta2, u = ul, v = vl, add = TRUE, by = 6)


mrf.exact.nll(par = theta, crf = fit, instances = samps)
mrf.nll(par = theta, crf = fit, instances = samps, infer.method = infer.exact)
?mrf.nll
fit$gradient


#------- USING CRF Library:
# Make a 2D plot of the gradient for the negative log likelihood:
uf <- CRF.grad.nll.mat.u
vf <- CRF.grad.nll.mat.v

# Sample negative log pseudo-likelihood:
image2D(x = theta1, y = theta2, z = CRF.nll.mat, contour = TRUE, xlab="theta1", ylab="theta2")
#image2D(x = theta1, y = theta2, z = nll.mat, contour = TRUE, xlab="theta1", ylab="theta2")


# Gradient:
quiver2D(x = theta1, y = theta2, u = uf, v = vf, add = TRUE, by = 6)


# Optimize theta: find minimum of negative log pseudo-likelihood:
# Auxiliary, gradient convenience function. Follows train.mrf in CRF:
gradient <- function(par, crf, ...) { crf$gradient }

fit$par <- c(0,0)                  # Reset theta to start at 0
opt.info  <- stats::optim(         # optimize parameters
  par          = fit$par,          # theta
  fn           = neglogpseudolik2, # objective function
  gr           = gradient,         # grad of obj func
  crf          = fit,              # passed to fn/gr
  samples      = samps,            # passed to fn/gr
  conditional.energy.function.type = "feature", # passed to fn/gr
  ff           = f0,               # passed to fn/gr
  gradQ        = TRUE,             # passed to fn/gr
  update.crfQ  = TRUE,             # passed to fn/gr
  method       = "L-BFGS-B",
  control      = list(trace = 1, REPORT=1))

opt.info$convergence
opt.info$message
fit$gradient
fit$nll
fit$par

# Make sure the potentials are generated from the optimized theta:
shift.pots(fit)
fit$node.pot
fit$edge.pot

# Instantiate a new CRF object to do the negative log likelihood minimization:
fit2 <- make.crf(adj, n.states)
fit2 <- make.features(fit2)
fit2 <- make.par(fit2, 2)
fit2$node.par[1,1,] <- 1
fit2$node.par[2,1,] <- 1

fit2$edge.par[[1]][1,1,1] <- 2
fit2$edge.par[[1]][2,2,1] <- 2
fit2$par.stat <- mrf.stat(crf = fit2, instances = samps)
fit2$par <- c(0,0)

opt.info  <- stats::optim(         # optimize parameters
  par          = fit2$par,         # theta
  fn           = negloglik,        # objective function
  gr           = gradient,         # grad of obj func
  crf          = fit2,             # passed to fn/gr
  samples      = samps,            # passed to fn/gr
  update.crfQ  = TRUE,             # passed to fn/gr
  method       = "L-BFGS-B",
  control      = list(trace = 1, REPORT=1))

opt.info$convergence
opt.info$message
fit2$gradient
fit2$nll
fit2$par
fit$par

fit2$par - fit$par

# Instantiate the true model with the actual theta:
tru <- make.crf(adj, n.states)
tru <- make.features(tru)
tru <- make.par(tru, 2)
tru$node.par[1,1,] <- 1
tru$node.par[2,1,] <- 1

tru$edge.par[[1]][1,1,1] <- 2
tru$edge.par[[1]][2,2,1] <- 2
tru$par <- c(1.098612, -0.7096765)

shift.pots(crf = tru)
shift.pots(crf = fit)
shift.pots(crf = fit2)

tru$node.pot
fit$node.pot
fit2$node.pot
tru$edge.pot
fit$edge.pot
fit2$edge.pot

tru.bel1  <- infer.exact(tru)
fit.bel1  <- infer.exact(fit)
fit2.bel1 <- infer.exact(fit2)

# Re-scale potentials to ensure they are in "standard form":
junk <- make.pots(parms = tru$par,  crf = tru,  rescaleQ = T, replaceQ = T)
junk <- make.pots(parms = fit$par,  crf = fit,  rescaleQ = T, replaceQ = T)
junk <- make.pots(parms = fit2$par, crf = fit2, rescaleQ = T, replaceQ = T)

tru$node.pot
fit$node.pot
fit2$node.pot
tru$edge.pot
fit$edge.pot
fit2$edge.pot

tru.bel2  <- infer.exact(tru)
fit.bel2  <- infer.exact(fit)
fit2.bel2 <- infer.exact(fit2)

tru.bel1$node.bel
tru.bel2$node.bel
tru.bel1$edge.bel
tru.bel2$edge.bel

fit.bel1$node.bel
fit.bel2$node.bel
fit.bel1$edge.bel
fit.bel2$edge.bel

fit2.bel1$node.bel
fit2.bel2$node.bel
fit2.bel1$edge.bel
fit2.bel2$edge.bel
