library(CRFutil)

grphf <- ~A:B
adj   <- ug(grphf, result="matrix")
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}

n.states <- 2
known.model <- make.crf(adj, n.states)

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


library(rgl)

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
n2p <- nodes2params.list2(fit, storeQ = T)


# Set up a grid of theta's overwhich to compute the sample pseudolikelihood:
theta1 <- seq(from =  2.0, to=2.5, by = 0.01)
theta2 <- seq(from =  0.4, to=0.8, by = 0.01)
theta.grid <- as.matrix(expand.grid(theta1,theta2))
dim(theta.grid)

# Compute the sample pseudolikelihoods over the theta grid:
neg.log.psls <- array(NA, nrow(theta.grid))
for(i in 1:nrow(theta.grid)) {
  print(i)
  theta           <- theta.grid[i,]
  neg.log.psls[i] <- neglogpseudolik(par = theta, crf = fit, samples = samps, ff = f0, update.crfQ = TRUE)

}

# Look for the rough minimum:
hist(neg.log.psls)
min(neg.log.psls)
min.idx <- which(neg.log.psls == min(neg.log.psls))
# Rough minimum:
c(theta.grid[min.idx,1], theta.grid[min.idx,2], neg.log.psls[min.idx])

# Plot the sample pseudo-likelihood over the theta grid and the rough minimum:
plot3d(theta.grid[,1], theta.grid[,2], neg.log.psls, xlab="theta1", ylab="theta2", zlab="")
text3d(theta.grid[min.idx,1], theta.grid[min.idx,2], neg.log.psls[min.idx], texts = "min")
lines3d(c(theta.grid[min.idx,1],theta.grid[min.idx,1]), c(theta.grid[min.idx,2],theta.grid[min.idx,2]), c(neg.log.psls[min.idx], max(neg.log.psls)))
text3d(theta.grid[min.idx,1], theta.grid[min.idx,2], max(neg.log.psls), texts = "min")


# Gradient of the sample neg. log pseudo-likelihood:
# Compute and plot gradient:
npll.mat        <- array(NA, c(length(theta1), length(theta2)))
grad.npll.mat.u <- array(NA, c(length(theta1), length(theta2)))
grad.npll.mat.v <- array(NA, c(length(theta1), length(theta2)))

count <- 1
for(i in 1:length(theta1)) {
  for(j in 1:length(theta2)) {

    print(count)
    theta            <- c(theta1[i], theta2[j])
    npll.mat[i,j]    <- neglogpseudolik(par = theta, crf = fit, samples = samps, ff = f0, update.crfQ = TRUE)
    neg.log.psl.grad <- grad.neglogpseudolik(par = theta, crf = fit, samples = samps, ff = f0)

    grad.npll.mat.u[i,j] <- neg.log.psl.grad[1]
    grad.npll.mat.v[i,j] <- neg.log.psl.grad[2]

    count <- count + 1
  }
}

# Make a 2D plot of the gradient:
library(OceanView)
u <- grad.npll.mat.u
v <- grad.npll.mat.v

# Sample negative log pseudo-likelihood:
image2D(x = theta1, y = theta2, z = npll.mat, contour = TRUE, xlab="theta1", ylab="theta2")

# Gradient:
quiver2D(x = theta1, y = theta2, u = u, v = v, add = TRUE, by = c(8,6))



# Nail down the min:
min(neg.log.psls)
min.idx <- which(neg.log.psls == min(neg.log.psls))

epsl <- 0.01
min.idx
theta1 <- seq(from = theta.grid[min.idx,1] - epsl, to=theta.grid[min.idx,1] + epsl, by = 0.0001)
theta2 <- seq(from = theta.grid[min.idx,2] - epsl, to=theta.grid[min.idx,2] + epsl, by = 0.0001)

npll.mat        <- array(NA, c(length(theta1), length(theta2)))
grad.npll.mat.u <- array(NA, c(length(theta1), length(theta2)))
grad.npll.mat.v <- array(NA, c(length(theta1), length(theta2)))
theta.mat       <- array(NA, c(length(theta1), length(theta2)))
prod(dim(theta.mat))

count <- 1
for(i in 1:length(theta1)) {
  for(j in 1:length(theta2)) {

    print(count)
    theta.mat[i,j]   <- count
    theta            <- c(theta1[i], theta2[j])
    npll.mat[i,j]    <- neglogpseudolik(par = theta, crf = fit, samples = samps, ff = f0, update.crfQ = TRUE)
    neg.log.psl.grad <- grad.neglogpseudolik(par = theta, crf = fit, samples = samps, ff = f0)

    grad.npll.mat.u[i,j] <- neg.log.psl.grad[1]
    grad.npll.mat.v[i,j] <- neg.log.psl.grad[2]

    count <- count + 1
  }
}

theta.grid2 <- as.matrix(expand.grid(theta1,theta2))
neg.log.psls2 <- as.numeric(npll.mat)
grad.neg.log.psls2.u <- as.numeric(grad.npll.mat.u)
grad.neg.log.psls2.v <- as.numeric(grad.npll.mat.v)

min(neg.log.psls2)
min.idx2 <- which(neg.log.psls2 == min(neg.log.psls2))
# Tighter minimum:
c(theta.grid2[min.idx2,1], theta.grid2[min.idx2,2], neg.log.psls2[min.idx2])
#Gradient at tighter min:
c(grad.neg.log.psls2.u[min.idx2], grad.neg.log.psls2.v[min.idx2])

# Plot the sample pseudo-likelihood over the theta grid and the rough minimum:
plot3d(theta.grid2[,1], theta.grid2[,2], neg.log.psls2, xlab="theta1", ylab="theta2", zlab="")
text3d(theta.grid2[min.idx2,1], theta.grid2[min.idx2,2], neg.log.psls2[min.idx2], texts = "min")
lines3d(c(theta.grid2[min.idx2,1],theta.grid2[min.idx2,1]), c(theta.grid2[min.idx2,2],theta.grid2[min.idx2,2]), c(neg.log.psls2[min.idx2], max(neg.log.psls2)))
text3d(theta.grid2[min.idx2,1], theta.grid2[min.idx2,2], max(neg.log.psls2), texts = "min")

# Color 3D sample negative log pseudo-likelihood sirface
nbcol <- 256
color <- rev(rainbow(nbcol, start = 0/6, end = 2/6)) #Color band width
zcol  <- cut(npll.mat, nbcol)
persp3d(theta1, theta2, npll.mat, aspect=c(1,1,1), col=color[zcol], xlab="",ylab="",zlab="")
lines3d(c(theta.grid2[min.idx2,1],theta.grid2[min.idx2,1]), c(theta.grid2[min.idx2,2],theta.grid2[min.idx2,2]), c(neg.log.psls2[min.idx2], max(neg.log.psls2)))
text3d(theta.grid2[min.idx2,1], theta.grid2[min.idx2,2], max(neg.log.psls2), texts = "min")


#
theta.opt <- c(theta.grid2[min.idx2,1], theta.grid2[min.idx2,2])
fit.pots <- make.pots(parms = theta.opt, crf = fit, rescaleQ=T, replaceQ=FALSE, printQ=FALSE)


known.model <- make.features(known.model)
known.model <- make.par(known.model, 2)
known.model$node.par[1,1,] <- 1
known.model$node.par[2,1,] <- 1
known.model$node.par

known.model$edge.par[[1]][1,1,1] <- 2
known.model$edge.par[[1]][2,2,1] <- 2
known.model$edge.par

known.theta <- make.par.from.potentials(known.model)
known.theta
known.pots <- make.pots(parms = known.theta, crf = known.model, rescaleQ=T, replaceQ=FALSE, printQ=FALSE)

# Ehhhh on the node potentials:
fit.pots[[1]]
known.pots[[1]]

# ????
fit.pots[[2]]
known.pots[[2]]

fit$par <- theta.opt
fit.pots <- make.pots(parms = theta.opt, crf = fit, rescaleQ=F, replaceQ=T, printQ=FALSE)
fit$par
fit$edge.pot
fit$node.pot

infer.exact(fit)
infer.exact(known.model)
# Not so great. Lets see what the MLE comes out to be...

# MLE fit:
fit2 <- make.crf(adj, n.states)
fit2 <- make.features(fit2)
fit2 <- make.par(fit2, 2)
fit2$node.par[1,1,] <- 1
fit2$node.par[2,1,] <- 1
fit2$node.par

fit2$edge.par[[1]][1,1,1] <- 2
fit2$edge.par[[1]][2,2,1] <- 2
fit2$edge.par

# First compute the sufficient stats needed by the likelihood and its’ grad
fit2$par.stat <- mrf.stat(fit2, samps)
fit2$par.stat

# Auxiliary, gradient convenience function. Follows train.mrf in CRF:
gradient <- function(par, crf, ...) { crf$gradient }

infr.meth <- infer.exact        # inference method needed for Z and marginals calcs
opt.info  <- stats::optim(    # optimize parameters
  par          = fit2$par,       # theta
  fn           = negloglik,     # objective function
  gr           = gradient,      # grad obj func
  crf          = fit2,           # passed to fn/gr
  samples      = samps,         # passed to fn/gr
  infer.method = infr.meth,     # passed to fn/gr
  update.crfQ  = TRUE,          # passed to fn/gr
  method       = "L-BFGS-B",
  control      = list(trace = 1, REPORT=1))
opt.info$convergence
opt.info$message
fit2$gradient
fit2$nll

infer.exact(fit2)
fit2$node.pot
fit2$edge.pot

fit2$par
fit$par

fit2$nll
min(neg.log.psls2)

fit2$nll
min(neg.log.psls)

fit2$par
c(theta.grid[min.idx,1], theta.grid[min.idx,2])
