library(CRFutil)

# Slayer field:
grphf <- ~A:B+A:C+A:D+A:E+B:C+B:D+B:E+C:D+D:E
gp    <- ug(grphf, result = "graph")
dev.off()
iplot(gp)

adj   <- ug(grphf, result="matrix")
adj

n.states     <- 2
slay        <- sim.field.random(adjacentcy.matrix=adj, num.states=n.states, num.sims=1000, seed=1)
samps       <- slay$samples
known.model <- slay$model


fit <- make.crf(adj, n.states)
fit <- make.features(fit)
fit <- make.par(fit, 14)

# One param per node:
for(i in 1:nrow(fit$node.par)){
  fit$node.par[i,1,1] <- i
}
fit$node.par

# Edge parameterization:
fit$edges # Check edge order first!
fit$edge.par[[1]][1,1,1] <- 6
fit$edge.par[[1]][2,2,1] <- 6
fit$edge.par[[2]][1,1,1] <- 7
fit$edge.par[[2]][2,2,1] <- 7
fit$edge.par[[3]][1,1,1] <- 8
fit$edge.par[[3]][2,2,1] <- 8
fit$edge.par[[4]][1,1,1] <- 9
fit$edge.par[[4]][2,2,1] <- 9
fit$edge.par[[5]][1,1,1] <- 10
fit$edge.par[[5]][2,2,1] <- 10
fit$edge.par[[6]][1,1,1] <- 11
fit$edge.par[[6]][2,2,1] <- 11
fit$edge.par[[7]][1,1,1] <- 12
fit$edge.par[[7]][2,2,1] <- 12
fit$edge.par[[8]][1,1,1] <- 13
fit$edge.par[[8]][2,2,1] <- 13
fit$edge.par[[9]][1,1,1] <- 14
fit$edge.par[[9]][2,2,1] <- 14

fit$edge.par

# Fit model to samples from the known model and obtain an estimate for w:
train.mrf(fit, samps, nll=mrf.exact.nll, infer.method = infer.exact)
fit$par.stat
w <- fit$par
w

# Shift potentials to w:
node.pot.orig <- fit$node.pot
edge.pot.orig <- fit$edge.pot
shift.pots(fit)
fit$node.pot
node.pot.orig

fit$edge.pot
edge.pot.orig
w

# Define states and feature function:
s1 <- 1
s2 <- 2
f  <- function(y){ as.numeric(c((y==s1),(y==s2))) }

# Enumerate all the state configurations
all.configs <- expand.grid(c(s1,s2),c(s1,s2),c(s1,s2),c(s1,s2),c(s1,s2))
colnames(all.configs) <- gp@nodes
all.configs


# Gen phi's for all possible configs:
M.all  <- compute.model.matrix(all.configs, fit$edges, fit$node.par, fit$edge.par, f)
M.all
colSums(M.all)

# Check sufficient stats for sample
M.samps  <- compute.model.matrix(samps, fit$edges, fit$node.par, fit$edge.par, f)
M.samps
colSums(M.samps) # Check
fit$par.stat     # Check

# All configuration energies:
fit.params.info   <- make.gRbase.potentials(fit, node.names=gp@nodes, state.nmes=c(s1,s2))
fit.params.info
w
fit.node.en       <- fit.params.info$node.energies
fit.edge.en       <- fit.params.info$edge.energies
fit.node.pot      <- fit.params.info$node.potentials
fit.edge.pot      <- fit.params.info$edge.potentials

dist.en.info <- distribution.from.energies(
  state.space = all.configs,
  edges.mat = fit$edges,
  node.energies = fit.node.en,
  edge.energies = fit.edge.en,
  energy.func = config.energy,
  ff = f)
dist.en.info


# All configuration energies computed with energy function and feature functions:
# Define a convenience function wrapper:
ener.func <- function(config) {
  engy <-config.energy(config = config, edges.mat = fit$edges,
                       one.lgp = fit.node.en,
                       two.lgp = fit.edge.en,
                       ff = f)
  return(engy)
}

config.energies <- sapply(1:nrow(all.configs), function(xx){ener.func(all.configs[xx,])})
prodPots        <- exp(config.energies)
Z               <- sum(prodPots)
Prs             <-prodPots/Z
Prs                          # Check
dist.en.info$state.probs     # Check
Prs-dist.en.info$state.probs # Check

# Compare config energies with energies computed as theta \dot phi:
alt.config.energies <- sapply(1:nrow(M.all), function(xx){w%*%M.all[xx,]})
config.energies                       # Check
alt.config.energies                   # Check
config.energies - alt.config.energies # Check

cbind(config.energies, alt.config.energies, config.energies-alt.config.energies)

