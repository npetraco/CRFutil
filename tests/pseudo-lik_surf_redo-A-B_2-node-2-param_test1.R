library(CRFutil)
library(rgl)

grphf <- ~A:B
adj   <- ug(grphf, result="matrix")
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}

n.states <- 2
known.model <- make.crf(adj, n.states)

# First lets get a sample of configurations from a "true" model:
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


neglogpseudolik(par = theta.grid[1,], crf = fit, samples = samps, ff = f0, update.crfQ = TRUE)

# Compute the sample pseudolikelihoods over the theta grid:
neg.log.psls <- array(NA, nrow(theta.grid))
for(i in 1:nrow(theta.grid)) {
  print(i)
  theta.i           <- theta.grid[i,]
  neg.log.psls[i] <- neglogpseudolik(par = theta.i, crf = fit, samples = samps, ff = f0, update.crfQ = TRUE)

  #negloglik(par, crf, samples, infer.method = infer.exact, update.crfQ = TRUE)

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
quiver2D(x = theta1, y = theta2, u = u, v = v, add = TRUE, by = 6)
