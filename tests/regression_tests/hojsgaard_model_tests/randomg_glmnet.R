library(igraph)
library(CRFutil)
library(glmnet)
library(rstan)
library(shinystan)
library(MASS)

# Make up a random graph
num.nodes <- 16
# START WITH FULLY SATURATED (pairwise) MODEL:
g         <- erdos.renyi.game(num.nodes, 1, typ="gnp")
dev.off()
plot(g)

# Get its adjacency matrix and genrate an MRF sample
adj <- as.matrix(as_adj(g))

f0       <- function(y){ as.numeric(c((y==1),(y==2)))}
rmod     <- make.empty.field(adj.mat = adj, parameterization.typ = "standard", plotQ = F)

# "true" theta
rmod$par <- runif(rmod$n.par,-1.5,1.5)
rmod$par

edgepar2edge <- cbind((num.nodes+1):rmod$n.par, rmod$edges)
colnames(edgepar2edge) <- c("par#","edgei","edgej")
edgepar2edge

# Knock out a few edges:
ko.num  <- 57 # knock out this many edges
ko.idxs <- sample((num.nodes+1):rmod$n.par, size = ko.num, replace = F)
rmod$par[ko.idxs] <- 0
rmod$par

# Make true pots from true theta
out.pot <- make.pots(parms = rmod$par,  crf = rmod,  rescaleQ = T, replaceQ = T)

# So now sample from the model as if we obtained an experimental sample:
num.samps <- 500
samps <- sample.exact(rmod, num.samps)
colnames(samps) <- 1:ncol(samps)
mrf.sample.plot(samps)

# Get frequency-counts of the configuration states:
configs.and.counts <- as.data.frame(ftable(data.frame(samps)))
freq.idx           <- ncol(configs.and.counts)
freqs              <- configs.and.counts[,freq.idx]
configs            <- configs.and.counts[,-freq.idx]

# Construct MRF model matrix.
# First get formula for graph
gf <- adj2formula(adj)
dev.off()
plot(ug(gf))

# Instantiate an empty model
glmn.fit   <- make.empty.field(graph.eq = gf, parameterization.typ = "standard")

print("Building model matrix")
MM  <- compute.model.matrix(
  configs   = configs,
  edges.mat = glmn.fit$edges,
  node.par  = glmn.fit$node.par,
  edge.par  = glmn.fit$edge.par,
  ff        = f0)
print("Done with model matrix. Sorry it's slow...")
dim(MM)
length(glmn.fit$par)

# Add intercept column
#MM <- cbind(rep(1,nrow(MM)), MM)


fit = glmnet(MM, freqs, family = "poisson")
print(fit)
plot(fit, label=T)

# These are the knocked out edges and their parameter numbers:
cbind(edgepar2edge[ko.idxs-num.nodes,], rmod$par[ko.idxs])

fit2 <- cv.glmnet(MM, freqs, family = "poisson")
plot(fit2)
fit2$lambda.min
log(fit2$lambda.min)

fit2$lambda.1se
log(fit2$lambda.1se)


# There are truely this many non zero theta:
rmod$n.par - ko.num
rmod$par

indic.true <- rep(1,rmod$n.par)
indic.true[ko.idxs] <- 0
indic.true <- c(1,indic.true)
length(indic.true)
indic.true


coef(fit2, s="lambda.min")
dim(coef(fit2, s="lambda.min"))

rmod.inf <- infer.exact(rmod)
rmod.inf$logZ
alp      <- log(nrow(configs)) - rmod.inf$logZ

cbind(coef(fit2, s="lambda.min"), indic.true, c(alp,rmod$par))
cbind(coef(fit2, s="lambda.1se"), indic.true, c(alp,rmod$par))

dim(MM)



