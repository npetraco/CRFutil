library(CRFutil)

# Put together MRF model:
grphf <- ~A:B + B:C + C:A
gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)

adj <- ug(grphf, result="matrix")
adj

n.states     <- 2
triangle.mrf <- make.crf(adj, n.states)
triangle.mrf <- make.features(triangle.mrf)
triangle.mrf <- make.par(triangle.mrf, 6)

triangle.mrf$node.par[1,1,1]      <- 1 # Represents tauA
triangle.mrf$node.par[2,1,1]      <- 2 # Represents tauB
triangle.mrf$node.par[3,1,1]      <- 3 # Represents tauC

triangle.mrf$edges # Check edge order first!
triangle.mrf$edge.par[[1]][1,1,1] <- 4 # Represents omegaAB
triangle.mrf$edge.par[[1]][2,2,1] <- 4
triangle.mrf$edge.par[[2]][1,1,1] <- 5 # Represents omegaAC
triangle.mrf$edge.par[[2]][2,2,1] <- 5
triangle.mrf$edge.par[[3]][1,1,1] <- 6 # Represents omegaBC
triangle.mrf$edge.par[[3]][2,2,1] <- 6

triangle.mrf$node.par
triangle.mrf$edge.par

# Define states and feature function:
s1 <- 1
s2 <- 2
f  <- function(y){ as.numeric(c((y==s1),(y==s2))) }


# Enumerate all the state configurations
config.mat <- expand.grid(c(s1,s2),c(s1,s2),c(s1,s2))
colnames(config.mat) <- c("A","B","C")
config.mat

# phi features:
phi.features(config = c(1,1,1), triangle.mrf$edges, triangle.mrf$node.par, triangle.mrf$edge.par, f)
phi.features(config = c(1,2,1), triangle.mrf$edges, triangle.mrf$node.par, triangle.mrf$edge.par, f)

phi.features.explicit(config = c(1,1,1), triangle.mrf$edges, triangle.mrf$node.par, triangle.mrf$edge.par, f)
phi.features.explicit(config = c(1,2,1), triangle.mrf$edges, triangle.mrf$node.par, triangle.mrf$edge.par, f)

# Model matrix and s:
cfgs <- rbind(c(1,1,1),c(1,2,1),c(1,1,1))
M    <- compute.model.matrix(cfgs, triangle.mrf$edges, triangle.mrf$node.par, triangle.mrf$edge.par, f)
s.M  <- colSums(M)
s.M

# Lets get some w
# First make a triangle model with known parameters:
known.model <- make.crf(adj, n.states)

# node weights:
tauA <- c( 1,    -1.3) + 1.3
tauA
tauB <- c(-0.85, -2.4) + 2.4
tauB
tauC <- c(3.82,   1.4) - 1.4
tauC
PsiA <- exp(tauA)
PsiB <- exp(tauB)
PsiC <- exp(tauC)

known.model$node.pot <- rbind(
  PsiA,
  PsiB,
  PsiC
)
known.model$node.pot

# edge weights:
omegaAB <- rbind(
  c( 3.5, -1.4) + 1.4,
  c(-1.4,  3.5) + 1.4
)
omegaAB
PsiAB <- exp(omegaAB)

omegaBC <- rbind(
  c( 2.6, 0.4) - 0.4,
  c( 0.4, 2.6) - 0.4
)
omegaBC
PsiBC <- exp(omegaBC)
PsiBC

omegaAC <- rbind(
  c(-0.6,  1.2) - 1.2,
  c( 1.2, -0.6) - 1.2
)
omegaAC
PsiAC <- exp(omegaAC)

known.model$edges # Check!
known.model$edge.pot[[1]] <- PsiAB
known.model$edge.pot[[2]] <- PsiAC
known.model$edge.pot[[3]] <- PsiBC
known.model$edge.pot

# Sample from the known model
samps <- sample.exact(known.model, 10000)

# Fit triangle model to samples from the known model and obtain an estimate for w:
train.mrf(triangle.mrf, samps, nll=mrf.exact.nll, infer.method = infer.exact)
triangle.mrf$par.stat
w <- triangle.mrf$par

triangle.mrf$node.pot
known.model$node.pot  # Huh??

triangle.mrf$edge.pot
known.model$edge.pot


w
log(triangle.mrf$node.pot)
w
log(triangle.mrf$edge.pot[[1]])
w
log(triangle.mrf$edge.pot[[2]])
w
log(triangle.mrf$edge.pot[[3]])


# Shift potentials to w:
triangle.mrf$node.pot
triangle.mrf$edge.pot
shift.pots(triangle.mrf)
triangle.mrf$node.pot
triangle.mrf$edge.pot


# Gen phi's for all possible configs:
all.configs <- expand.grid(c(1,2),c(1,2),c(1,2))
M.all  <- compute.model.matrix(all.configs, triangle.mrf$edges, triangle.mrf$node.par, triangle.mrf$edge.par, f)
M.all

# All configuration energies:
triangle.params.info   <- make.gRbase.potentials(triangle.mrf, node.names=gp@nodes, state.nmes=c(s1,s2))
triangle.node.en  <- triangle.params.info$node.energies
triangle.edge.en  <- triangle.params.info$edge.energies
triangle.node.pot <- triangle.params.info$node.potentials
triangle.edge.pot <- triangle.params.info$edge.potentials

dist.en.info <- distribution.from.energies(
  state.space = all.configs,
  edges.mat = triangle.mrf$edges,
  node.energies = triangle.node.en,
  edge.energies = triangle.edge.en,
  energy.func = config.energy,
  ff = f)
dist.en.info


# All configuration energies computed with energy function and feature functions:
# Define a convenience function wrapper:
ener.func <- function(config) {
  engy <-config.energy(config = config, edges.mat = triangle.mrf$edges,
                       one.lgp = triangle.node.en,
                       two.lgp = triangle.edge.en,
                       ff = f)
  return(engy)
}

config.energies <- sapply(1:nrow(all.configs), function(xx){ener.func(all.configs[xx,])})
prodPots        <- exp(config.energies)
Z               <- sum(prodPots)
Prs             <-prodPots/Z
Prs                       # Check
dist.en.info$state.probs  # Check

# Compare config energies with energies computed as theta \dot phi:
alt.config.energies <- sapply(1:nrow(M.all), function(xx){w%*%M.all[xx,]})
config.energies                       # Check
alt.config.energies                   # Check
config.energies - alt.config.energies # Check

