# Fit triangle model all ways with obscene number of obs
# Generate exact and empirical data for comparison with fitting methods here:

library(CRFutil)

# Path to save data to:
fpth <- "/home/npetraco/codes/R/CRFutil/tests/regression_tests/hojsgaard_model_tests/triangle_model/triangle_data/"

# Triangle model
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
num.samps <- 3500
set.seed(1)
samps <- sample.exact(true.model, num.samps)
mrf.sample.plot(samps)

# ******* Uses the names to watch for node name switching amongs other fitting methods ******
node.names      <- c("X.1","X.2","X.3")
colnames(samps) <- node.names

# Write sample to file for use with all fitting methods:
save(samps, file = paste0(fpth,"triangle_samp.RData"))
# To load the data again
#load(paste0(fpth,"triangle_samp.RData"))

# Generate and save the true state probabilities
# ******** Use these node names instead or gp@nodes: **********
gp@nodes
node.names

pot.info.true.model        <- make.gRbase.potentials(true.model, node.names = node.names, state.nmes = c("1","2"))
gR.dist.info.true.model    <- distribution.from.potentials(pot.info.true.model$node.potentials,
                                                           pot.info.true.model$edge.potentials)
#logZ.true.model            <- gR.dist.info.true.model$logZ
joint.dist.info.true.model <- as.data.frame(as.table(gR.dist.info.true.model$state.probs))
joint.dist.info.true.model
# Rearrange node order. It's probably out of wack
exact.dist.info <- joint.dist.info.true.model[,c(2,3,1,4)]
colnames(exact.dist.info)[4] <- "Exact.Freq"
exact.dist.info
save(samps, file = paste0(fpth,"triangle_exact_dist.RData"))

# Fold samples into contingency table for a look
X.cont <- xtabs(~., data=data.frame(samps))
X.cont

# Take a look at the empirical frequency-counts of the states:
X.counts <- as.data.frame(ftable(data.frame(samps)))
sum(X.counts[,4]) # Should be the sample size
num.samps == sum(X.counts[,4])

# Estimate state probabilities is by the EMPIRICAL relative frquencies:
X.freq <- X.counts
X.freq[,4] <- X.freq[,4]/sum(X.counts[,4])
X.freq <- cbind(X.freq, X.counts[,4]) # Carry along the counts the computation is based on
colnames(X.freq)[c(4,5)] <- c("Emp.Freq","Counts")
X.freq
save(X.freq, file = paste0(fpth,"triangle_empirical_dist.RData"))
