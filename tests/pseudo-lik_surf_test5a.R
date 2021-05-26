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
theta1 <- seq(from =   0.0, to=3.5, by = 0.05)
theta2 <- seq(from =  -2.0, to=2.0, by = 0.05)
theta.grid <- as.matrix(expand.grid(theta1,theta2))
dim(theta.grid)

# Compute the sample pseudolikelihoods over the theta grid:
#neg.log.psls <- array(NA, nrow(theta.grid))
neg.log.lik <- array(NA, nrow(theta.grid))
fit$par.stat <- mrf.stat(crf = fit, instances = samps)
for(i in 1:nrow(theta.grid)) {
  print(i)
  theta           <- theta.grid[i,]
  #neg.log.psls[i] <- neglogpseudolik(par = theta, crf = fit, samples = samps, ff = f0, update.crfQ = TRUE)
  neg.log.lik[i] <- mrf.nll(par = theta, crf = fit, instances = samps, infer.method = infer.exact)
}

# Look for the rough minimum:
hist(neg.log.lik)
min(neg.log.lik)
min.idx.lik <- which(neg.log.lik == min(neg.log.lik))
# Rough minimum:
c(theta.grid[min.idx.lik,1], theta.grid[min.idx.lik,2], neg.log.lik[min.idx.lik])

# Plot the sample pseudo-likelihood over the theta grid and the rough minimum:
plot3d(theta.grid[,1], theta.grid[,2], neg.log.lik, xlab="theta1", ylab="theta2", zlab="")
text3d(theta.grid[min.idx.lik,1], theta.grid[min.idx.lik,2], neg.log.lik[min.idx.lik], texts = "min")
lines3d(c(theta.grid[min.idx.lik,1],theta.grid[min.idx.lik,1]), c(theta.grid[min.idx.lik,2],theta.grid[min.idx.lik,2]), c(neg.log.lik[min.idx.lik], max(neg.log.lik)))
text3d(theta.grid[min.idx.lik,1], theta.grid[min.idx.lik,2], max(neg.log.lik), texts = "min")
