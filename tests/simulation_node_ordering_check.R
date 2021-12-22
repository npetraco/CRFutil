# Check ordering of nodes then simulations are run!:

library(CRFutil)

# Graph:
grphf1 <- ~X.1:X.2 + X.2:X.3 + X.3:X.4 + X.4:X.1
grphf2 <- ~X1:X2 + X2:X3 + X3:X4 + X4:X1
grphf3 <- ~1:2 + 2:3 + 3:4 + 4:1
grphf4 <- ~1:4 + 4:2 + 2:3 + 3:1

# Type 1 node labels:
adj1                <- ug(grphf1, result="matrix")
node2nme1           <- data.frame(1:ncol(adj1),colnames(adj1))
colnames(node2nme1) <- c("node","name")
node.names1         <- colnames(adj1)
node.names1
gp1 <- ug(grphf1, result = "graph")
dev.off()
iplot(gp1)


# Type 2 node labels:
adj2                <- ug(grphf2, result="matrix")
node2nme2           <- data.frame(1:ncol(adj2),colnames(adj2))
colnames(node2nme2) <- c("node","name")
node.names2         <- colnames(adj2)
node.names2
gp2 <- ug(grphf2, result = "graph")
dev.off()
iplot(gp2)

# Type 3 node labels:
adj3                <- ug(grphf3, result="matrix")
node2nme3           <- data.frame(1:ncol(adj3),colnames(adj3))
colnames(node2nme3) <- c("node","name")
node.names3         <- colnames(adj3)
node.names3
gp3 <- ug(grphf3, result = "graph")
dev.off()
iplot(gp3)

# Type 4 node labels:
adj4                <- ug(grphf4, result="matrix")
node2nme4           <- data.frame(1:ncol(adj4),colnames(adj4))
colnames(node2nme4) <- c("node","name")
node.names4         <- colnames(adj4)
node.names4
gp4 <- ug(grphf4, result = "graph")
dev.off()
iplot(gp4)

#
model1.info <- sim.field.random(adjacentcy.matrix=adj1, num.states=2, num.sims=10000, seed=3102)
samps1      <- model1.info$samples
samps1

model2.info <- sim.field.random(adjacentcy.matrix=adj2, num.states=2, num.sims=10000, seed=3102)
samps2      <- model2.info$samples
samps2

model3.info <- sim.field.random(adjacentcy.matrix=adj3, num.states=2, num.sims=10000, seed=3102)
samps3      <- model3.info$samples
samps3

model4.info <- sim.field.random(adjacentcy.matrix=adj4, num.states=2, num.sims=10000, seed=3102)
samps4      <- model4.info$samples
samps4

SA <- samps1
SB <- samps3
sum(!sapply(1:nrow(SA), function(xx){!sum(!(SA[xx,] == SB[xx,]))}))

# So, I'd doesn't matter which way we label the nodes as long as the conectivity is consistent.

# But haw can we test what the actual column names of the sample matrix are?? Do the nodes get permuted at some point??
# I guess we can take a hudge sample and see if it can reproduce the node and edge beliefs.

model1a.info <- sim.field.random(adjacentcy.matrix=adj1, num.states=2, num.sims=1000000, seed=3102)
samps1a <- model1a.info$samples

fit <- make.empty.field(adj.mat = adj1, parameterization.typ = "standard", plotQ = T)

# MLE for parameters of model. Follows train.mrf in CRF, just gives more control and output:
gradient      <- function(par, crf, ...) { crf$gradient } # Auxiliary, gradient convenience function.
fit$par.stat <- mrf.stat(fit, samps1a)   # requisite sufficient statistics
infr.meth <- infer.exact               # inference method needed for Z and marginals calcs
opt.info  <- stats::optim(             # optimize parameters
  par          = fit$par,              # theta
  fn           = negloglik,            # objective function
  gr           = gradient,             # grad of obj func
  crf          = fit,                  # passed to fn/gr
  samples      = samps1a,                # passed to fn/gr
  infer.method = infr.meth,            # passed to fn/gr
  update.crfQ  = TRUE,                 # passed to fn/gr
  method       = "L-BFGS-B",
  control      = list(trace = 1, REPORT=1))

# Checks: May have to re-run optimize a few times to get gradient down:
opt.info$convergence
opt.info$message
fit$gradient
fit$nll
fit$par         # Estimated parameter vector

bels.true <- infer.junction(model1a.info$model)
bels.fit  <- infer.junction(fit)

bels.fit$node.bel
bels.true$node.bel
bels.fit$edge.bel[[1]]
bels.true$edge.bel[[1]]

bels.fit$edge.bel[[2]]
bels.true$edge.bel[[2]]

bels.fit$edge.bel[[3]]
bels.true$edge.bel[[3]]

bels.fit$edge.bel[[4]]
bels.true$edge.bel[[4]]

# Beliefs are approx well by fit to sample so I guess CRF isn't changing to order of the nodes compatred to adj
