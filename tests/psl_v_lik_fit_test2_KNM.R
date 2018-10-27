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

knm <- make.crf(adj, n.states)
knm <- make.features(knm)
knm <- make.par(knm, 15)
knm$node.par[1,1,] <- 1
knm$node.par[2,1,] <- 2
knm$node.par[3,1,] <- 3
knm$node.par[4,1,] <- 4
knm$node.par[5,1,] <- 5
knm$node.par

knm$edges
knm$edge.par[[1]][1,1,1] <- 6
knm$edge.par[[1]][2,2,1] <- 6
knm$edge.par[[2]][1,1,1] <- 7
knm$edge.par[[2]][2,2,1] <- 7
knm$edge.par[[3]][1,1,1] <- 8
knm$edge.par[[3]][2,2,1] <- 8
knm$edge.par[[4]][1,1,1] <- 9
knm$edge.par[[4]][2,2,1] <- 9
knm$edge.par[[5]][1,1,1] <- 10
knm$edge.par[[5]][2,2,1] <- 10
knm$edge.par[[6]][1,1,1] <- 11
knm$edge.par[[6]][2,2,1] <- 11
knm$edge.par[[7]][1,1,1] <- 12
knm$edge.par[[7]][2,2,1] <- 12
knm$edge.par[[8]][1,1,1] <- 13
knm$edge.par[[8]][2,2,1] <- 13
knm$edge.par[[9]][1,1,1] <- 14
knm$edge.par[[9]][2,2,1] <- 14
knm$edge.par[[10]][1,1,1] <- 15
knm$edge.par[[10]][2,2,1] <- 15
knm$edge.par

#knm$par <- runif(15,-1.5,1.1)
#knm$par # "true" theta
knm$par <- c(-0.15932366, 0.32110660, 0.31562139, -0.08118330, 0.52367998, 0.57853082, -0.82745675, -0.05998841, -0.63907303, 0.98135824, -0.56745243, -0.18572940, -0.03587448, 0.93756418, 0.43754708)
out.pot <- make.pots(parms = knm$par,  crf = knm,  rescaleQ = T, replaceQ = T)
knm$node.pot
knm$edge.pot

# So now sample from the model as if we obtained an experimental sample:
num.samps <- 25
#set.seed(1)
samps <- sample.exact(knm, num.samps)
mrf.sample.plot(samps)
write.csv(samps, file = "/Users/npetraco/latex/papers/dust/steph_diss/CRFutil/tests/psl_v_lik_fit_test2_KNM.csv")

knm.bel <- infer.exact(knm)
knm.bel$node.bel
# [,1]      [,2]
# [1,] 0.3898024 0.6101976
# [2,] 0.5956328 0.4043672
# [3,] 0.6679497 0.3320503
# [4,] 0.4907574 0.5092426
# [5,] 0.6697415 0.3302585

knm.bel$edge.bel
# [[1]]
# [,1]      [,2]
# [1,] 0.2629105 0.1268919
# [2,] 0.3327223 0.2774753
#
# [[2]]
# [,1]      [,2]
# [1,] 0.1748250 0.2149774
# [2,] 0.4931247 0.1170729
#
# [[3]]
# [,1]      [,2]
# [1,] 0.1593869 0.2304155
# [2,] 0.3313705 0.2788271
#
# [[4]]
# [,1]      [,2]
# [1,] 0.1663788 0.2234236
# [2,] 0.5033627 0.1068349
#
# [[5]]
# [,1]      [,2]
# [1,] 0.4575297 0.1381031
# [2,] 0.2104200 0.1939472
#
# [[6]]
# [,1]      [,2]
# [1,] 0.2251281 0.3705047
# [2,] 0.2656294 0.1387379
#
# [[7]]
# [,1]      [,2]
# [1,] 0.3909200 0.2047128
# [2,] 0.2788215 0.1255457
#
# [[8]]
# [,1]      [,2]
# [1,] 0.3300835 0.3378662
# [2,] 0.1606739 0.1713764
#
# [[9]]
# [,1]      [,2]
# [1,] 0.5447001 0.1232496
# [2,] 0.1250414 0.2070089
#
# [[10]]
# [,1]      [,2]
# [1,] 0.3754370 0.1153204
# [2,] 0.2943045 0.2149381

# Distribution calculations:
s1 <- 1
s2 <- 2
# Enumerate all the state configurations
config.mat <- expand.grid(c(s1,s2),c(s1,s2),c(s1,s2),c(s1,s2),c(s1,s2))
colnames(config.mat) <- c("1","2","3","4","5")
config.mat

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
  ff            = f0)
knm.Pr.en <- knm.dist.en.info$state.probs
knm.Pr.en
# Known model configuration probabilities for testing:
c(0.039406231, 0.119251532, 0.012757382, 0.122790885, 0.010001263, 0.005784078, 0.023048807, 0.042396813,
  0.053563908, 0.143769266, 0.005574245, 0.047586668, 0.012653248, 0.006490468, 0.009373736, 0.015292983,
  0.013462865, 0.011348665, 0.003006157, 0.008059803, 0.022283599, 0.003589821, 0.035420642, 0.018148862,
  0.043902920, 0.032824296, 0.003151263, 0.007493629, 0.067636492, 0.009664142, 0.034559644, 0.015705688)


sum(knm.Pr.en)
cbind(config.mat,round(100*knm.Pr.en))

# As a check compute Probs of all configurations and logZ from potentials as well:
knm.dist.pot.info <- distribution.from.potentials(
  gRbase.node.potentials = knm.node.pot,
  gRbase.edge.potentials = knm.edge.pot)
knm.Pr.pot <- knm.dist.pot.info$state.probs
knm.Pr.pot
sum(knm.Pr.pot)
