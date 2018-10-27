library(CRFutil)

# Fully connected
grphf <- ~1:2+1:3+1:4+1:5+2:3+2:4+2:5+3:4+3:5+4:5
adj <- ug(grphf, result="matrix")
# Check the graph:
gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)

f0 <- function(y){ as.numeric(c((y==1),(y==2)))}
n.states <- 2

# Read in configuration data:
samps <- read.csv(file = "/Users/npetraco/latex/papers/dust/steph_diss/CRFutil/tests/psl_v_lik_fit_test2_KNM.csv")
samps

# Instantiate an empty model to fit:
lik <- make.crf(adj, n.states)
lik <- make.features(lik)
lik <- make.par(lik, 15)
lik$node.par[1,1,] <- 1
lik$node.par[2,1,] <- 2
lik$node.par[3,1,] <- 3
lik$node.par[4,1,] <- 4
lik$node.par[5,1,] <- 5
lik$edge.par[[1]][1,1,1] <- 6
lik$edge.par[[1]][2,2,1] <- 6
lik$edge.par[[2]][1,1,1] <- 7
lik$edge.par[[2]][2,2,1] <- 7
lik$edge.par[[3]][1,1,1] <- 8
lik$edge.par[[3]][2,2,1] <- 8
lik$edge.par[[4]][1,1,1] <- 9
lik$edge.par[[4]][2,2,1] <- 9
lik$edge.par[[5]][1,1,1] <- 10
lik$edge.par[[5]][2,2,1] <- 10
lik$edge.par[[6]][1,1,1] <- 11
lik$edge.par[[6]][2,2,1] <- 11
lik$edge.par[[7]][1,1,1] <- 12
lik$edge.par[[7]][2,2,1] <- 12
lik$edge.par[[8]][1,1,1] <- 13
lik$edge.par[[8]][2,2,1] <- 13
lik$edge.par[[9]][1,1,1] <- 14
lik$edge.par[[9]][2,2,1] <- 14
lik$edge.par[[10]][1,1,1] <- 15
lik$edge.par[[10]][2,2,1] <- 15


# Optimize theta: find minimum of negative log pseudo-likelihood:
# Auxiliary, gradient convenience function. Follows train.mrf in CRF:
gradient <- function(par, crf, ...) { crf$gradient }

# likelihood fit:
lik$par.stat <- mrf.stat(crf = lik, instances = as.matrix(samps))
#M.mat <- compute.model.matrix(samps, lik$edges, lik$node.par, lik$edge.par, f0)
#lik$par.stat <- colSums(M.mat)
lik$par
lik$par.stat


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

true.theta <- c(-0.15932366, 0.32110660, 0.31562139, -0.08118330, 0.52367998, 0.57853082, -0.82745675, -0.05998841, -0.63907303, 0.98135824, -0.56745243, -0.18572940, -0.03587448, 0.93756418, 0.43754708)
cbind(lik$par,true.theta)


# Ensure stored potentials conform to fit theta
lik$node.pot
lik$edge.pot
junk.lik <- make.pots(parms = lik$par,  crf = lik,  rescaleQ = T, replaceQ = T)
lik$node.pot
lik$edge.pot

lik.bel <- infer.exact(lik)

lik.bel$node.bel
lik.bel$edge.bel

# Distribution calculations:
s1 <- 1
s2 <- 2
# Enumerate all the state configurations
config.mat <- expand.grid(c(s1,s2),c(s1,s2),c(s1,s2),c(s1,s2),c(s1,s2))
colnames(config.mat) <- c("1","2","3","4","5")
config.mat

# Full distribution likelihood estimate of theta
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

# Likelihood distribution from pseudolikelihood estimate of theta
psl.Pr.en <- c(
  0.0022575398, 0.0396155182, 0.0030344525, 0.1071287100, 0.0013680352, 0.0060647305, 0.0003165662, 0.0028234123, 0.1377301222, 0.1487460208,
  0.0384408167, 0.0835227201, 0.0522888870, 0.0142662595, 0.0025124379, 0.0013790870, 0.0011473678, 0.0058576769, 0.0076216177, 0.0782826526,
  0.0100562542, 0.0129700997, 0.0115001373, 0.0298404724, 0.0182861934, 0.0057455628, 0.0252224091, 0.0159437851, 0.1004096614, 0.0079701972,
  0.0238430068, 0.0038075901)

# Pseudolikelihood distribution
psl2.Pr.en <- c(
  1.531268e-04, 4.058400e-02, 2.979772e-04, 2.245701e-01, 1.721896e-04, 7.025222e-03, 5.365273e-06, 4.259462e-04, 2.366274e-01, 2.306907e-01,
  3.612705e-02, 8.911863e-02, 6.883615e-02, 7.700862e-03, 1.537353e-04, 4.427305e-05, 4.366533e-05, 1.177818e-03, 5.101413e-03, 2.153183e-01,
  1.465426e-02, 4.961460e-02, 2.826395e-02, 1.125004e-01, 5.433729e-03, 4.885701e-04, 5.555158e-02, 6.233371e-03, 3.785510e-01, 3.945161e-03,
  4.907314e-02, 7.129123e-04)


# Known model configuration probabilities for testing:
true.config.Prs <- c(
  0.039406231, 0.119251532, 0.012757382, 0.122790885, 0.010001263, 0.005784078, 0.023048807, 0.042396813,
  0.053563908, 0.143769266, 0.005574245, 0.047586668, 0.012653248, 0.006490468, 0.009373736, 0.015292983,
  0.013462865, 0.011348665, 0.003006157, 0.008059803, 0.022283599, 0.003589821, 0.035420642, 0.018148862,
  0.043902920, 0.032824296, 0.003151263, 0.007493629, 0.067636492, 0.009664142, 0.034559644, 0.015705688)

trup    <- round(100*true.config.Prs)
likp    <- round(100*lik.Pr.en)
psllikp <- round(100*psl.Pr.en)
pslp    <- round(100*psl2.Pr.en)
cbind(config.mat,trup,likp,psllikp,pslp)


# Examine differences between the estimated distributions and the true distribution:
zero.eps <- runif(1,0,1) * 1e-5
trup.mod <- trup
trup.mod[which(trup==0)] <- zero.eps
trup.mod

likp.mod <- likp
likp.mod[which(likp==0)] <- zero.eps
likp.mod

psllikp.mod <- psllikp
psllikp.mod[which(psllikp==0)] <- zero.eps
psllikp.mod

pslp.mod <- pslp
pslp.mod[which(pslp==0)] <- zero.eps
pslp.mod

# KL divergences:
trup.mod/likp.mod
trup.mod * log(trup.mod/likp.mod)
likp.kld <- sum(trup.mod * log(trup.mod/likp.mod))

trup.mod/psllikp.mod
trup.mod * log(trup.mod/psllikp.mod)
psllikp.kld <- sum(trup.mod * log(trup.mod/psllikp.mod))

trup.mod/pslp.mod
trup.mod * log(trup.mod/pslp.mod)
pslp.kld <- sum(trup.mod * log(trup.mod/pslp.mod))

likp.kld
psllikp.kld
pslp.kld
