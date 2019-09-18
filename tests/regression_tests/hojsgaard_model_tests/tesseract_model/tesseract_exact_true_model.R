#library(rstan)
library(CRFutil)

# Tesseract field model:
grphf <-
  ~X.1:X.2   + X.1:X.4   + X.1:X.5  + X.1:X.13 +
  X.2:X.3   + X.2:X.6   + X.2:X.14 +
  X.3:X.4   + X.3:X.7   + X.3:X.15 +
  X.4:X.8   + X.4:X.16  +
  X.5:X.6   + X.5:X.8   + X.5:X.9 +
  X.6:X.7   + X.6:X.10  +
  X.7:X.8   + X.7:X.11  +
  X.8:X.12  +
  X.9:X.10  + X.9:X.12  + X.9:X.13 +
  X.10:X.11 + X.10:X.14 +
  X.11:X.12 + X.11:X.15 +
  X.12:X.16 +
  X.13:X.14 + X.13:X.16 +
  X.14:X.15 +
  X.15:X.16

adj <- ug(grphf, result="matrix")
adj
node2nme <- data.frame(1:ncol(adj),colnames(adj))
colnames(node2nme) <- c("node","name")
node2nme

# Load reference sample data
#setwd("/home/npetraco/codes/R/CRFutil/tests/regression_tests/hojsgaard_model_tests/tesseract_model/")
setwd("tests/regression_tests/hojsgaard_model_tests/tesseract_model/")
load(file = "samps.RData")

# Regenerate the "true" model
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}
tess <- make.empty.field(adj.mat = adj, parameterization.typ = "standard", plotQ = T)
set.seed(87)
tess$par <- runif(tess$n.par,-1.5,1.5)
tess$par # "true" theta
out.pot <- make.pots(parms = tess$par,  crf = tess,  rescaleQ = F, replaceQ = T)
tess$edges
tess$node.pot
tess$edge.pot

node.names      <- colnames(adj)
node.names   # Check
node2nme[,2] # Check
colnames(samps) <- node.names
head(samps)
mrf.sample.plot(samps)


gR.tess     <- make.gRbase.potentials(tess, node.names = node.names, state.nmes = c("1","2"))

gR.dist.info.fit.model <- distribution.from.potentials(gR.tess$node.potentials, gR.tess$edge.potentials)
logZ.fit.model         <- gR.dist.info.fit.model$logZ
logZ.fit.model # log partition function from BP with a theta (from regression)

pri <- gR.dist.info.fit.model$state.probs
dim(pri)
joint.dist.info.fit.model <- as.data.frame(as.table(gR.dist.info.fit.model$state.probs))
dim(joint.dist.info.fit.model)

# Reaggange column states to be in the same order as node.names
scrambled.nn <- colnames(joint.dist.info.fit.model)[1:(ncol(joint.dist.info.fit.model)-1)]
scrambled.nn
node.names

rearr.idxs <- sapply(1:length(node.names), function(xx){which(scrambled.nn == node.names[xx])})
scrambled.nn[rearr.idxs]
node.names

joint.dist.info.fit.model.rearr <- joint.dist.info.fit.model[c(rearr.idxs,ncol(joint.dist.info.fit.model))]
colnames(joint.dist.info.fit.model.rearr)

joint.configs <- joint.dist.info.fit.model.rearr[,-ncol(joint.dist.info.fit.model.rearr)]
joint.configs <- as.matrix(joint.configs)
joint.configs[which(joint.configs == 2, arr.ind = T)] <- 0
head(joint.configs) # Columns of the joint dist should now be in node.names order a

# Compute index of each configuration in the node.names order and BINARY (character)
nn <- length(node.names)
powers.of.two <- 2^(0:(nn - 1))

joint.config.idxs <- sapply(1:nrow(joint.configs), function(xx){as.integer(as.numeric(joint.configs[xx,]) %*% powers.of.two)})
dim(joint.configs)
length(joint.config.idxs)

# Order joint.config.idxs into decreasing order
# start 1 1 1 1 ...
# ends  0 0 0 0 ...
joint.configs.row.order <- order(joint.config.idxs, decreasing = T)
joint.configs.row.order

head(data.frame(joint.config.idxs[joint.configs.row.order], joint.configs[joint.configs.row.order,])) # Check

config.prs <- joint.dist.info.fit.model.rearr[,ncol(joint.dist.info.fit.model.rearr)]
plot(config.prs, typ="h")
plot(config.prs[joint.configs.row.order], typ="h", ylim=c(0,0.006)) # Re-order

joint.dist.info.ordered <- data.frame(joint.config.idxs[joint.configs.row.order], joint.configs[joint.configs.row.order,], config.prs[joint.configs.row.order])
colnames(joint.dist.info.ordered) <- c("config.idx", node.names,"config.Pr")
head(joint.dist.info.ordered) # Check

sum(config.prs[joint.configs.row.order]) # config probs sum to 1?

save(joint.dist.info.ordered, file = "tesseract_true_model.dist.RData")
