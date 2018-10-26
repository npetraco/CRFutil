library(CRFutil)

grphf <- ~1:2+1:3+2:3
adj <- ug(grphf, result="matrix")
# Check the graph:
gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)

f0 <- function(y){ as.numeric(c((y==1),(y==2)))}
n.states <- 2

knm <- make.crf(adj, n.states)
knm <- make.features(knm)
knm <- make.par(knm, 9)
knm$node.par[1,1,] <- 1
knm$node.par[2,1,] <- 2
knm$node.par[3,1,] <- 3
knm$node.par

knm$edges
knm$edge.par[[1]][1,1,1] <- 4
knm$edge.par[[1]][2,2,1] <- 5
knm$edge.par[[2]][1,1,1] <- 6
knm$edge.par[[2]][2,2,1] <- 7
knm$edge.par[[3]][1,1,1] <- 8
knm$edge.par[[3]][2,2,1] <- 9
knm$edge.par

knm$par <- runif(9,-3,3)
knm$par # "true" theta
# 0.7271790 -1.0614092 -1.2386324  0.8606001  1.6952350 -1.9150751 -1.5859374 -1.5238856 -1.8766993 # just incase
out.pot <- make.pots(parms = knm$par,  crf = knm,  rescaleQ = T, replaceQ = T)
#knm$node.pot
#knm$edge.pot

# So now sample from the model as if we obtained an experimental sample:
num.samps <- 25
#set.seed(1)
samps <- sample.exact(knm, num.samps)
mrf.sample.plot(samps)

# Instantiate a models to fit:
# For the pseudolikelihood
psl <- make.crf(adj, n.states)
psl <- make.features(psl)
psl <- make.par(psl, 9)
psl$node.par[1,1,] <- 1
psl$node.par[2,1,] <- 2
psl$node.par[3,1,] <- 3
psl$edge.par[[1]][1,1,1] <- 4
psl$edge.par[[1]][2,2,1] <- 5
psl$edge.par[[2]][1,1,1] <- 6
psl$edge.par[[2]][2,2,1] <- 7
psl$edge.par[[3]][1,1,1] <- 8
psl$edge.par[[3]][2,2,1] <- 9
psl$par

# For the likelihood
lik <- make.crf(adj, n.states)
lik <- make.features(lik)
lik <- make.par(lik, 9)
lik$node.par[1,1,] <- 1
lik$node.par[2,1,] <- 2
lik$node.par[3,1,] <- 3
lik$edge.par[[1]][1,1,1] <- 4
lik$edge.par[[1]][2,2,1] <- 5
lik$edge.par[[2]][1,1,1] <- 6
lik$edge.par[[2]][2,2,1] <- 7
lik$edge.par[[3]][1,1,1] <- 8
lik$edge.par[[3]][2,2,1] <- 9
lik$par.stat <- mrf.stat(crf = lik, instances = samps)
lik$par
lik$par.stat


# Optimize theta: find minimum of negative log pseudo-likelihood:
# Auxiliary, gradient convenience function. Follows train.mrf in CRF:
gradient <- function(par, crf, ...) { crf$gradient }

# likelihood fit:
opt.info  <- stats::optim(         # optimize parameters
  par          = lik$par,          # theta
  fn           = negloglik,        # objective function
  gr           = gradient,         # grad of obj func
  crf          = lik,             # passed to fn/gr
  samples      = samps,            # passed to fn/gr
  update.crfQ  = TRUE,             # passed to fn/gr
  method       = "L-BFGS-B",
  control      = list(trace = 1, REPORT=1))

opt.info$convergence
opt.info$message
lik$gradient
lik$nll
lik$par
knm$par

# pseudolikelihood fit:
opt.info  <- stats::optim(         # optimize parameters
  par          = psl$par,          # theta
  fn           = neglogpseudolik2, # objective function
  gr           = gradient,         # grad of obj func
  crf          = psl,              # passed to fn/gr
  samples      = samps,            # passed to fn/gr
  conditional.energy.function.type = "feature", # passed to fn/gr
  ff           = f0,               # passed to fn/gr
  gradQ        = TRUE,             # passed to fn/gr
  update.crfQ  = TRUE,             # passed to fn/gr
  method       = "L-BFGS-B",
  control      = list(trace = 1, REPORT=1))

opt.info$convergence
opt.info$message
psl$gradient
psl$nll
psl$par
knm$par
lik$par

# Make sure potentials conform with theta of their model
out.pot <- make.pots(parms = knm$par,  crf = knm,  rescaleQ = T, replaceQ = T)
out.pot <- make.pots(parms = psl$par,  crf = psl,  rescaleQ = T, replaceQ = T)
out.pot <- make.pots(parms = lik$par,  crf = lik,  rescaleQ = T, replaceQ = T)

knm.bel <- infer.exact(knm)
psl.bel <- infer.exact(psl)
lik.bel <- infer.exact(lik)

knm.bel$node.bel
psl.bel$node.bel
lik.bel$node.bel

knm.bel$edge.bel
psl.bel$edge.bel
lik.bel$edge.bel



# Distribution calculations:
s1 <- 1
s2 <- 2
# Enumerate all the state configurations
config.mat <- expand.grid(c(s1,s2),c(s1,s2),c(s1,s2))
colnames(config.mat) <- c("1","2","3")




# Full distribution but pseudo likelihood estimate of theta
# Get energies from potentials anf decorate both with gRBase annotations:
psl.params   <- make.gRbase.potentials(crf=psl, node.names=gp@nodes, state.nmes=c(s1,s2))
psl.node.en  <- psl.params$node.energies
psl.edge.en  <- psl.params$edge.energies
psl.node.pot <- psl.params$node.potentials
psl.edge.pot <- psl.params$edge.potentials

# Compute Probs of all configurations and logZ from energies:
psl.dist.en.info <- distribution.from.energies(
  state.space   = config.mat,
  edges.mat     = psl$edges,
  node.energies = psl.node.en,
  edge.energies = psl.edge.en,
  energy.func   = config.energy,
  ff = f0)
psl.Pr.en <- psl.dist.en.info$state.probs
psl.Pr.en
cbind(config.mat,round(100*psl.Pr.en))

# As a check compute Probs of all configurations and logZ from potentials as well:
psl.dist.pot.info <- distribution.from.potentials(
  gRbase.node.potentials = psl.node.pot,
  gRbase.edge.potentials = psl.edge.pot)
psl.Pr.pot <- psl.dist.pot.info$state.probs
psl.Pr.pot






# Full distribution with respect to "true" theta
# Get energies from potentials anf decorate both with gRBase annotations:
knm.params   <- make.gRbase.potentials(crf=knm, node.names=gp@nodes, state.nmes=c(s1,s2))
knm.node.en  <- knm.params$node.energies
knm.edge.en  <- knm.params$edge.energies
knm.node.pot <- knm.params$node.potentials
knm.edge.pot <- knm.params$edge.potentials

# Compute Probs of all configurations and logZ from energies:
knm.dist.en.info <- distribution.from.energies(
  state.space   = config.mat,
  edges.mat     = knm$edges,
  node.energies = knm.node.en,
  edge.energies = knm.edge.en,
  energy.func   = config.energy,
  ff = f0)
knm.Pr.en <- knm.dist.en.info$state.probs
knm.Pr.en
cbind(config.mat,round(100*knm.Pr.en), round(100*psl.Pr.en))

# As a check compute Probs of all configurations and logZ from potentials as well:
knm.dist.pot.info <- distribution.from.potentials(
  gRbase.node.potentials = knm.node.pot,
  gRbase.edge.potentials = knm.edge.pot)
knm.Pr.pot <- knm.dist.pot.info$state.probs
knm.Pr.pot





# Full distribution estimate of theta
# Get energies from potentials anf decorate both with gRBase annotations:
lik.params   <- make.gRbase.potentials(crf=lik, node.names=gp@nodes, state.nmes=c(s1,s2))
lik.node.en  <- lik.params$node.energies
lik.edge.en  <- lik.params$edge.energies
lik.node.pot <- lik.params$node.potentials
lik.edge.pot <- lik.params$edge.potentials

# Compute Probs of all configurations and logZ from energies:
lik.dist.en.info <- distribution.from.energies(
  state.space   = config.mat,
  edges.mat     = lik$edges,
  node.energies = lik.node.en,
  edge.energies = lik.edge.en,
  energy.func   = config.energy,
  ff = f0)
lik.Pr.en <- lik.dist.en.info$state.probs
lik.Pr.en
cbind(config.mat, round(100*knm.Pr.en), round(100*lik.Pr.en), round(100*psl.Pr.en))

# As a check compute Probs of all configurations and logZ from potentials as well:
lik.dist.pot.info <- distribution.from.potentials(
  gRbase.node.potentials = lik.node.pot,
  gRbase.edge.potentials = lik.edge.pot)
lik.Pr.pot <- lik.dist.pot.info$state.probs
lik.Pr.pot






# Full pseudo likelihood distribution estimate using psl optimized theta
psl2.dist.pot.info <- pseudolikelihoods.from.energies2(
  state.space  = config.mat,
  crf          = psl,
  cond.en.form = "feature",
  ff           = f0)
# psl2.dist.pot.info$condtional.energies
# psl2.dist.pot.info$complement.condtional.energies
# psl2.dist.pot.info$conditional.Zs
# psl2.dist.pot.info$conditional.Prs
# psl2.dist.pot.info$complement.conditional.Prs

psl2.Pr.en <- psl2.dist.pot.info$pseudo.likelihoods
psl2.Pr.en
cbind(config.mat, round(100*knm.Pr.en), round(100*lik.Pr.en), round(100*psl.Pr.en), round(100*psl2.Pr.en))
