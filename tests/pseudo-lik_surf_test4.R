library(CRFutil)
library(rgl)
library(OceanView)

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
    c(3.2, 6.1),
    c(4.3, 3.2)
  )

known.model$node.pot[1,]  <- PsiA
known.model$node.pot[2,]  <- PsiB
known.model$edge.pot[[1]] <- PsiAB

# So now sample from the model as if we obtained an experimental sample:
num.samps <- 25
set.seed(4)
samps <- sample.exact(known.model, num.samps)
mrf.sample.plot(samps)

#----------------------------------------------------------------------
# Now fit a model using a sample from the known model.
# For now we determine theta by hand. We can do this by looking at the
# negative log pseudolikelihood surface because there are only two parameters:

fit2 <- make.crf(adj, n.states)
fit2 <- make.features(fit2)
fit2 <- make.par(fit2, 2)
#fit2$node.par[1,1,] <- 1
#fit2$node.par[2,1,] <- 1
fit2$node.par

fit2$edge.par[[1]][1,1,1] <- 1
fit2$edge.par[[1]][1,2,1] <- 2
fit2$edge.par[[1]][2,2,1] <- 1
fit2$edge.par
n2p <- nodes2params.list2(fit2, storeQ = T) # Generalized to account to nodes without params...????
n2p

theta1 <- seq(from =  -1.5, to=0.75, by = 0.05)
theta2 <- seq(from =  -1.5, to=1.0, by = 0.05)

epsl <- 0.01
theta1 <- seq(from = theta.grid[min.idx,1] - epsl, to=theta.grid[min.idx,1] + epsl, by = 0.0003)
theta2 <- seq(from = theta.grid[min.idx,2] - epsl, to=theta.grid[min.idx,2] + epsl, by = 0.0003)


length(theta1)
length(theta2)
length(theta1)*length(theta2)
theta.grid <- as.matrix(expand.grid(theta1,theta2))
dim(theta.grid)

# Now lets try looping over the different parameter vectors contained in theta.grid. For each compute the
# corresponding sample pseudolikelihood:
neg.log.psls      <- array(NA, nrow(theta.grid))
grad.neg.log.psls <- array(NA, c(nrow(theta.grid), 2))
for(i in 1:nrow(theta.grid)) {

  print(i)
  theta                 <- theta.grid[i,]
  neg.log.psls[i]       <- neglogpseudolik(par = theta, crf = fit2, samples = samps, ff = f0, update.crfQ = TRUE)
  grad.neg.log.psls[i,] <- grad.neglogpseudolik(par = theta, crf = fit2, samples = samps, ff = f0)

}
neg.log.psls
hist(neg.log.psls)
min(neg.log.psls)


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

min.idx <- which(neg.log.psls == min(neg.log.psls))
min.idx
c(theta.grid[min.idx,1], theta.grid[min.idx,2], neg.log.psls[min.idx])
grad.neg.log.psls[min.idx,] # Gradient at the min


plot3d(theta.grid[,1], theta.grid[,2], neg.log.psls)
text3d(theta.grid[min.idx,1], theta.grid[min.idx,2], neg.log.psls[min.idx], texts = "min")
lines3d(c(theta.grid[min.idx,1],theta.grid[min.idx,1]), c(theta.grid[min.idx,2],theta.grid[min.idx,2]), c(neg.log.psls[min.idx], max(neg.log.psls)))


# Plot gradient:
npll.mat2        <- array(NA, c(length(theta1), length(theta2)))
grad.npll.mat.u2 <- array(NA, c(length(theta1), length(theta2)))
grad.npll.mat.v2 <- array(NA, c(length(theta1), length(theta2)))

count <- 1
for(i in 1:length(theta1)) {
  for(j in 1:length(theta2)) {

    print(count)
    theta <- c(theta1[i], theta2[j])

    npll.mat2[i,j] <- neglogpseudolik(par = theta, crf = fit2, samples = samps, ff = f0, update.crfQ = TRUE)
    neg.log.psl.grad2 <- grad.neglogpseudolik(par = theta, crf = fit2, samples = samps, ff = f0)

    grad.npll.mat.u2[i,j] <- neg.log.psl.grad2[1]
    grad.npll.mat.v2[i,j] <- neg.log.psl.grad2[2]

    count <- count + 1
  }
}

image2D(x = theta1, y = theta2, z = npll.mat2, contour = TRUE)
u2 <- grad.npll.mat.u2
v2 <- grad.npll.mat.v2
quiver2D(x = theta1, y = theta2, u = u2, v = v2, add = TRUE, by = 3)

# min theta, nlps-llk:
c(theta.grid[min.idx,1], theta.grid[min.idx,2], neg.log.psls[min.idx])
