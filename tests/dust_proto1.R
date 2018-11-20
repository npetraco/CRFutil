library(CRFutil)
library(rstan)
library(shinystan)
library(coda)

# Model:
grphf <- ~1:2 + 1:3 + 1:4 + 1:6 + 2:3 + 2:4 + 3:4 + 3:5 + 5:6 + 5:7
gamt <-  rbind(
  c(1,2),
  c(1,3),
  c(1,4),
  c(1,6),
  c(2,3),
  c(2,4),
  c(3,4),
  c(3,5),
  c(5,6),
  c(5,7))

adj <- edges2adj(gamt, plotQ = T)

knm <- make.empty.field(
  graph.eq             = NULL,
  adj.mat              = adj,
  parameterization.typ = "ising2",
  node.par             = NULL,
  edge.par             = NULL,
  plotQ                = T)
dump.crf(knm)
knm$edges


set.seed(6)
knm$par <- runif(6,-1.5,1.1)
knm$par # "true" theta
out.pot <- make.pots(parms = knm$par,  crf = knm,  rescaleQ = T, replaceQ = T)


# So now sample from the model as if we obtained an experimental sample:
num.samps <- 25
set.seed(1)
samps <- sample.exact(knm, num.samps)
mrf.sample.plot(samps)
