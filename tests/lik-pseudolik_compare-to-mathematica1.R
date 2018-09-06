library(CRFutil)
library(rgl)


#
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
samps
junk <- write.csv(samps, file = "/Users/npetraco/latex/papers/dust/steph_diss/CRFutil/tests/samps.csv")
#mrf.sample.plot(samps)

#----------------------------------------------------------------------
# Now fit a model using a sample from the known model.
# For now we determine theta by hand. We can do this by looking at the
# negative log pseudolikelihood surface because there are only two parameters:

fit2 <- make.crf(adj, n.states)
fit2 <- make.features(fit2)
fit2 <- make.par(fit2, 2)
fit2$node.par[1,1,] <- 1
fit2$node.par[2,1,] <- 1
fit2$node.par

fit2$edge.par[[1]][1,1,1] <- 2
fit2$edge.par[[1]][2,2,1] <- 2
fit2$edge.par
n2p <- nodes2params.list(fit2, storeQ = T)

theta1 <- seq(from =  2.0, to=2.5, by = 0.01)
theta2 <- seq(from =  0.4, to=0.8, by = 0.01)

#epsl <- 0.005
#theta1 <- seq(from = theta.grid[min.idx,1] - epsl, to=theta.grid[min.idx,1] + epsl, by = 0.0003)
#theta2 <- seq(from = theta.grid[min.idx,2] - epsl, to=theta.grid[min.idx,2] + epsl, by = 0.0003)

length(theta1)
length(theta2)
length(theta1)*length(theta2)
theta.grid <- as.matrix(expand.grid(theta1,theta2))
dim(theta.grid)

#-----------PSEUDO LOG LIKELIHOOD-------------------------
# Now lets try looping over the different parameter vectors contained in theta.grid. For each compute the
# corresponding sample pseudolikelihood:
num.theta         <- 5
neg.log.psls      <- array(NA, num.theta)
grad.neg.log.psls <- array(NA, c(num.theta, 2))
for(i in 1:num.theta) {

  print(i)
  theta                 <- theta.grid[i,]
  print(paste(theta[1], theta[2]))
  neg.log.psls[i]       <- neglogpseudolik(par = theta, crf = fit2, samples = samps, ff = f0, update.crfQ = FALSE)
  print(neg.log.psls[i])
  #grad.neg.log.psls[i,] <- grad.neglogpseudolik(par = theta, crf = fit2, samples = samps, ff = f0)

}
neg.log.psls
hist(neg.log.psls)
min(neg.log.psls)

theta.grid[1:5,]

#STOP HERE FOR NOW

#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

min.idx <- which(neg.log.psls == min(neg.log.psls))
min.idx
c(theta.grid[min.idx,1], theta.grid[min.idx,2], neg.log.psls[min.idx])
grad.neg.log.psls[min.idx,] # Gradient at the min


plot3d(theta.grid[,1], theta.grid[,2], neg.log.psls)
text3d(theta.grid[min.idx,1], theta.grid[min.idx,2], neg.log.psls[min.idx], texts = "min")
lines3d(c(theta.grid[min.idx,1],theta.grid[min.idx,1]), c(theta.grid[min.idx,2],theta.grid[min.idx,2]), c(neg.log.psls[min.idx], max(neg.log.psls)))


#-----------LOG LIKELIHOOD-------------------------
theta1 <- seq(from =  -10, to=10, by = 0.1)
theta2 <- seq(from =  -10, to=10, by = 0.1)

#epsl <- 0.005
#theta1 <- seq(from = theta.grid[min.idx,1] - epsl, to=theta.grid[min.idx,1] + epsl, by = 0.0003)
#theta2 <- seq(from = theta.grid[min.idx,2] - epsl, to=theta.grid[min.idx,2] + epsl, by = 0.0003)

length(theta1)
length(theta2)
length(theta1)*length(theta2)
theta.grid <- as.matrix(expand.grid(theta1,theta2))
dim(theta.grid)

neg.log.lik      <- array(NA, nrow(theta.grid))
grad.neg.log.lik <- array(NA, c(nrow(theta.grid), 2))
fit2$par.stat    <- mrf.stat(crf = fit2, instances = samps)
fit2$par.stat
for(i in 1:nrow(theta.grid)) {

  print(i)
  theta                 <- theta.grid[i,]
  neg.log.lik[i]        <- negloglik(par = theta, crf = fit2, samples = samps, infer.method = infer.exact, update.crfQ = F)
  #grad.neg.log.lik[i,] <- grad.negloglik(par = theta, crf = fit2, samples = samps, ff = f0)

}
neg.log.lik
hist(neg.log.lik)
min(neg.log.lik)
plot3d(theta.grid[,1], theta.grid[,2], neg.log.lik)
