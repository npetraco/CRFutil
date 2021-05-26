library(CRFutil)
#library(rstan)
#library(shinystan)
#library(coda)
#library(gRbase)
#library(gRim)
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
out.pot <- make.pots(parms = true.model$par,  crf = true.model,  rescaleQ = F, replaceQ = T)
infer.exact(true.model)

# So now sample from the model as if we obtained an experimental sample:
num.samps <- 25
set.seed(1)
samps <- sample.exact(true.model, num.samps)
mrf.sample.plot(samps)
colnames(samps) <- c("X.1","X.2","X.3")

# Raw case list (Observed states) to Aggregated Frequency Table
Xfreq <- as.data.frame(ftable(data.frame(samps)))
#ftable(data.frame(samps))

## Raw case-list (Observed states) to contingency table
Xcont <- xtabs(~., data=data.frame(samps))

# Contrasts to build Model matrix ??????
X1.dc <- Xfreq[,1]
contrasts(X1.dc)
contrasts(X1.dc) <- contr.sum(2)
contrasts(X1.dc) # Only really need this????

X2.dc <- Xfreq[,2]
contrasts(X2.dc)
contrasts(X2.dc) <- contr.sum(2)
contrasts(X2.dc) # Only really need this????

X3.dc <- Xfreq[,3]
contrasts(X3.dc)
contrasts(X3.dc) <- contr.sum(2)
contrasts(X3.dc) # Only really need this????

# Contrasts to build Model matrix ??????
Xm <- model.matrix(~X1.dc + X2.dc + X3.dc + X1.dc:X2.dc + X1.dc:X3.dc + X2.dc:X3.dc,
                   contrasts = list(X1.dc = "contr.sum",
                                    X2.dc = "contr.sum",
                                    X3.dc = "contr.sum"))
Xm

ll3 <- loglm(Freq~X1.dc + X2.dc + X3.dc + X1.dc:X2.dc + X1.dc:X3.dc + X2.dc:X3.dc, data = Xfreq)
coef(ll3)
Xfreq


ttm <- matrix(c(1,2),c(2,1))
rownames(ttm) <- c(1,2)
ttm
tt1 <- X1.dc
contrasts(tt1) <- ttm
tt1
tt2 <- X2.dc
contrasts(tt2) <- ttm
tt2
tt3 <- X3.dc
contrasts(tt3) <- ttm
tt3


ll4 <- loglm(Freq~tt1 + tt2 + tt3 + tt1:tt2 + tt1:tt3 + tt2:tt3, data = Xfreq)
ll4
# Look at emp relative freqs vs the fitted relative freqs:
X.fitted4 <- fitted(ll4)
X.Prob.fitted4 <- X.fitted4/sum(X.fitted4)
X.Prob.fitted4


# Look at emp relative freqs vs the fitted relative freqs:
X.fitted <- fitted(ll3)
X.Prob.fitted <- X.fitted/sum(X.fitted)
as.data.frame(ftable(X.Prob.fitted))
as.data.frame(as.table(X.Prob.fitted))

pot.info.true.model        <- make.gRbase.potentials(true.model, node.names = gp@nodes, state.nmes = c("1","2"))
gR.dist.info.true.model    <- distribution.from.potentials(pot.info.true.model$node.potentials, pot.info.true.model$edge.potentials)
logZ.true.model            <- gR.dist.info.true.model$logZ
joint.dist.info.true.model <- as.data.frame(as.table(gR.dist.info.true.model$state.probs))
ext <- joint.dist.info.true.model[,c(2,3,1,4)]

#
llm3 <- as.data.frame(as.table(X.Prob.fitted))
llm3

llm4 <- as.data.frame(as.table(X.Prob.fitted4))
llm4

sr.idxs <- sapply(1:nrow(ext), function(xx){row.match(ext[xx,1:3],llm3[,1:3])})
cbind(llm3[sr.idxs,],llm4[sr.idxs,], ext)

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

resul <- cbind(old[,1:3], llm3[sr.idxs,4]*100, llm4[sr.idxs,4]*100, old[,4:6], ext[,4]*100)
resul # ????????????????????

