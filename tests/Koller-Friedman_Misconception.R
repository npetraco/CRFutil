library(CRFutil)

# Put together the Koller-Freidman "Misconception" MRF model and get a sample from it:
grphf <- ~A:B + B:C + C:D + D:A
gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)

adj   <- ug(grphf, result="matrix")
adj

n.states <- 2
known.model <- make.crf(adj, n.states)

# Node pots not in original Misconception example: ????
PsiA <- c(1,1)
PsiB <- c(1,1)
PsiC <- c(1,1)
PsiD <- c(1,1)

# Edge pots in the Misconception example:
PsiAB <-
  rbind(
    c(30, 5),
    c(1, 10)
  )
PsiBC <-
  rbind(
    c(100, 1),
    c(1, 100)
  )
PsiCD <-
  rbind(
    c(1, 100),
    c(100, 1)
  )
PsiDA <-
  rbind(
    c(100, 1),
    c(1, 100)
  )

known.model$node.pot[1,] <- PsiA
known.model$node.pot[2,] <- PsiB
known.model$node.pot[3,] <- PsiC
known.model$node.pot[4,] <- PsiD
known.model$node.pot

known.model$edges # Check!
known.model$edge.pot[[1]] <- PsiAB
known.model$edge.pot[[2]] <- PsiDA
known.model$edge.pot[[3]] <- PsiBC
known.model$edge.pot[[4]] <- PsiCD

# Check again!
known.model$node.pot
known.model$edge.pot

# Check: Get same "Misconception" config probs?:
pot.info <- make.gRbase.potentials(known.model, node.names = gp@nodes, state.nmes = c("0","1"))
pot.info

gR.dist.info    <- distribution.from.potentials(pot.info$node.potentials, pot.info$edge.potentials)
logZ            <- gR.dist.info$logZ
joint.dist.info <- as.data.frame(as.table(gR.dist.info$state.probs))
joint.dist.info

joint.dist.info2 <- cbind(joint.dist.info[,c(3,4,2,1)], round( exp(logZ) * joint.dist.info[,5],6), round(joint.dist.info[,5],6))
colnames(joint.dist.info2) <- c("A","B","C","D", "Unnormalized", "Normalized")
joint.dist.info2 # Yup, we do.

# So now sample from the model as if we obtained an experimental sample:
num.samps <- 50
set.seed(1)
samps <- sample.exact(known.model, num.samps)

# Let's go banans and use the more flexible parameterization for the model we'll fit.
# Also we will assume the node parameters are not all 0. This will
# Require 1 param per node (4 node params) and 2 params per edge
# (2*4 = 8 edge params) for a total of 12 parameters.
# Initalize the fit model:
fit <- make.crf(adj, n.states)
fit <- make.features(fit)
fit <- make.par(fit, 12)
length(fit$par)

# Fill the parameter index matrices with the standard parameterization
# Standard parameterization, one param per node:
for(i in 1:nrow(fit$node.par)){
  fit$node.par[i,1,1] <- i
}
fit$node.par

# Use the more flexible edge parameterization:
fit$edges # Check edge order first!
fit$edge.par[[1]][1,1,1] <- 5
fit$edge.par[[1]][2,2,1] <- 6
fit$edge.par[[2]][1,1,1] <- 7
fit$edge.par[[2]][2,2,1] <- 8
fit$edge.par[[3]][1,1,1] <- 9
fit$edge.par[[3]][2,2,1] <- 10
fit$edge.par[[4]][1,1,1] <- 11
fit$edge.par[[4]][2,2,1] <- 12
fit$edge.par
#===================================


# Gradient function of the log-likelihood
grad.llk <- function(crf, inference.func, samples) {
  N        <- nrow(samples)
  dL       <- N*feature.means(crf, inference.func) - crf$par.stat
  #crf$grad <- dL
  return(dL)
}

# Test
fit$par.stat <- mrf.stat(fit, samps)
grad.llk(fit, infer.junction, samps)
# Same?
feature.means(fit, infer.junction)*num.samps - mrf.stat(fit, samps)

mrf.stat(fit, samps)  # First term of the gradient (its constant): N \times \hat{\text{E}}_{\boldsymbol \theta}[{\boldsymbol \phi}]
feature.means(fit, infer.junction)*num.samps


#--------------------
gradient <- function(par, fit, ...) { fit$gradient }
fit$par.stat <- mrf.stat(fit, samps)
opt.info <- stats::optim(fit$par, mrf.junction.nll, gradient, fit, samps, infer.junction,
             method = "L-BFGS-B",
             control = list(trace = 1, REPORT=1))
opt.info$convergence
opt.info$message

# Reset:
fit$par      <- numeric(length(fit$par))
fit$gradient <- numeric(length(fit$gradient))

mrf.junction.nll(fit$par, fit, samps, infer.junction)
fit$gradient
grad.llk(fit, infer.junction, samps)

# [1]  0.59534139 -0.04613655 -0.16954273  0.13159944  0.12603018 -0.42317465  1.96125385  1.23431303  3.40588695
# [10]  3.62156623 -2.57678713 -2.53884385
#mrf.update(fit)

#
# Check
train.mrf(fit, samps, nll=mrf.junction.nll, infer.method = infer.junction, trace=1)
fit$par
fit$gradient
grad.llk(fit, infer.junction, samps)
# [1]  0.59534139 -0.04613655 -0.16954273  0.13159944  0.12603018 -0.42317465  1.96125385  1.23431303  3.40588695
# [10]  3.62156623 -2.57678713 -2.53884385


#
#--------------------
gradient <- function(par, fit, ...) { fit$gradient }
fit$par.stat <- mrf.stat(fit, samps)
opt.info <- stats::optim(fit$par, mrf.exact.nll, gradient, fit, samps, infer.exact,
             method = "L-BFGS-B",
             control = list(trace = 1, REPORT=1))
opt.info$convergence
opt.info$message
# Reset:
fit$par      <- numeric(length(fit$par))
fit$gradient <- numeric(length(fit$gradient))

mrf.exact.nll(fit$par, fit, samps, infer.exact)
fit$gradient
grad.llk(fit, infer.lbp, samps)

#----------------------
opt.info <- stats::optim(fit$par, mrf.lbp.nll, gradient, fit, samps, infer.lbp,
                         method = "L-BFGS-B",
                         control = list(trace = 1, REPORT=1))
opt.info$convergence
opt.info$message
# Reset:
fit$par      <- numeric(length(fit$par))
fit$gradient <- numeric(length(fit$gradient))
