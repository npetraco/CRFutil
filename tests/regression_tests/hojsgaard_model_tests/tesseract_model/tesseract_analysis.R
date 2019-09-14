library(rstan)
library(shinystan)
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

gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)
plot(gp)

# Load essentials and results from Stan Poisson fit
setwd("/home/npetraco/codes/R/CRFutil/tests/regression_tests/hojsgaard_model_tests/tesseract_model/")
load(file = "bfit.RData")
load(file = "M.RData")
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
node.names
node2nme[,2]
colnames(samps) <- node.names
head(samps)
mrf.sample.plot(samps)

N     <- nrow(samps)
theta <- extract(bfit, "theta")[[1]]
alpha <- extract(bfit, "alpha")[[1]]

theta.est <- apply(X = theta, MARGIN = 2, FUN = median)
alpha.est <- median(alpha)
logZ.est  <- log(N) - alpha.est

# Configutation energies
E.est <- M %*% theta.est

# Configuration probabilities:
pr.est <- exp(-logZ.est + E.est)
plot(pr.est*100, typ="h")
sum(pr.est) # Missing some %

# Configuration indices of M?????

# "True model" configuration probabilities