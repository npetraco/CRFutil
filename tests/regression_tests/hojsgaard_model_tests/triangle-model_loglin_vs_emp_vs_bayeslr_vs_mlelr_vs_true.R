library(CRFutil)
#library(rstan)
#library(shinystan)
#library(coda)
#library(gRbase)
library(gRim)
library(MASS)



# Fully connected
grphf <- ~1:2+1:3+2:3
adj <- ug(grphf, result="matrix")
# Check the graph:
gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)

f0 <- function(y){ as.numeric(c((y==1),(y==2)))}
n.states <- 2

# Instantiate the "true model":
true.model <- make.crf(adj, n.states)
true.model <- make.features(true.model)
true.model <- make.par(true.model, 6)
true.model$node.par[1,1,] <- 1
true.model$node.par[2,1,] <- 2
true.model$node.par[3,1,] <- 3
true.model$edge.par[[1]][1,1,1] <- 4
true.model$edge.par[[1]][2,2,1] <- 4
true.model$edge.par[[2]][1,1,1] <- 5
true.model$edge.par[[2]][2,2,1] <- 5
true.model$edge.par[[3]][1,1,1] <- 6
true.model$edge.par[[3]][2,2,1] <- 6

set.seed(6)
true.model$par <- runif(6,-1.5,1.1)
true.model$par # "true" theta
out.pot <- make.pots(parms = true.model$par,  crf = true.model,  rescaleQ = T, replaceQ = T)


# So now sample from the model as if we obtained an experimental sample:
num.samps <- 25
set.seed(1)
samps <- sample.exact(true.model, num.samps)
mrf.sample.plot(samps)

# Fold the observed configurations into a 2^3 (#states^dim(configs) for Ising and Potts-like models) -way contingency table
colnames(samps) <- c("1","2","3")
X <- xtabs(~., data=samps)
X
dim(X)
ftable(X)
dimnames(X)

# Hojsgaard: One estimate of the config probabilities is by the relative frquencies:
X.Prob <- X/sum(X)
ftable(X.Prob)
as.data.frame(as.table(X.Prob))

# First fit the Hojgaard model with loglm and gRim
ll2 <- loglm(grphf, data=X); # loglm from MASS which uses loglin in base
ll2
X.loglm.coefs <- coef(ll2)
X.loglm.coefs

# Look at emp relative freqs vs the fitted relative freqs:
X.fitted <- fitted(ll2)
X.Prob.fitted <- X.fitted/sum(X.fitted)
ftable(X.Prob.fitted)
data.frame(as.data.frame(as.table(X.Prob)), as.data.frame(as.table(X.Prob.fitted)))

pot.info.true.model        <- make.gRbase.potentials(true.model, node.names = gp@nodes, state.nmes = c("1","2"))
gR.dist.info.true.model    <- distribution.from.potentials(pot.info.true.model$node.potentials, pot.info.true.model$edge.potentials)
logZ.true.model            <- gR.dist.info.true.model$logZ
joint.dist.info.true.model <- as.data.frame(as.table(gR.dist.info.true.model$state.probs))
joint.dist.info.true.model[,c(2,3,1,4)]

#
emp <- as.data.frame(as.table(X.Prob))
llm <- as.data.frame(as.table(X.Prob.fitted))
ext <- joint.dist.info.true.model[,c(2,3,1,4)]
#
sr.idxs <- sapply(1:nrow(ext), function(xx){row.match(ext[xx,1:3],emp[,1:3])})
cbind(emp[sr.idxs,], llm[sr.idxs,], ext)
#sapply(1:nrow(llm), function(xx){row.match(ext[xx,1:3],llm[,1:3])})

# From triangle_logistic-regress_bayes.R
old <- rbind(
  c(1,  1,  1,         10.3,       10.5,          18.4),
  c(1,  1,  2,          1.1,        1.5,           8.0),
  c(2,  1,  1,         22.0,       21.5,          15.6),
  c(2,  1,  2,         23.1,       22.5,          22.6),
  c(1,  2,  1,          9.1,        9.5,           4.2),
  c(1,  2,  2,         10.3,       10.5,          14.9),
  c(2,  2,  1,          2.0,        2.5,           1.3),
  c(2,  2,  2,         22.0,       21.5,          15.0)
)
colnames(old) <- c("1","2","3","bayes.lrm", "mle.lrm", "true.model") 

resul <- cbind(emp[sr.idxs,4]*100, llm[sr.idxs,4]*100, old[,4:6], ext[,4]*100)
resul
colnames(resul)[1:2] <- c("emp", "llm")
colnames(resul)[6] <- c("true.too")

resul
