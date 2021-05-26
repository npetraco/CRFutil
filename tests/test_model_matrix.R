library(CRFutil)

# Put together MRF model:
grphf <- ~A:B + B:C + C:A
gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)

adj <- ug(grphf, result="matrix")
adj

n.states     <- 2
triangle.mrf <- make.crf(adj, n.states)
triangle.mrf <- make.features(triangle.mrf)
triangle.mrf <- make.par(triangle.mrf, 6)

triangle.mrf$node.par[1,1,1]      <- 1 # Represents tauA
triangle.mrf$node.par[2,1,1]      <- 2 # Represents tauB
triangle.mrf$node.par[3,1,1]      <- 3 # Represents tauC

triangle.mrf$edges # Check edge order first!
triangle.mrf$edge.par[[1]][1,1,1] <- 4 # Represents omegaAB
triangle.mrf$edge.par[[1]][2,2,1] <- 4
triangle.mrf$edge.par[[2]][1,1,1] <- 5 # Represents omegaAC
triangle.mrf$edge.par[[2]][2,2,1] <- 5
triangle.mrf$edge.par[[3]][1,1,1] <- 6 # Represents omegaBC
triangle.mrf$edge.par[[3]][2,2,1] <- 6

triangle.mrf$node.par
triangle.mrf$edge.par

# Define states and feature function:
s1 <- 1
s2 <- 2
f  <- function(y){ as.numeric(c((y==s1),(y==s2))) }


# Enumerate all the state configurations
config.mat <- expand.grid(c(s1,s2),c(s1,s2),c(s1,s2))
colnames(config.mat) <- c("A","B","C")
config.mat

triangle.mrf$node.par
triangle.mrf$edge.par
phi.features(config = c(2,2,2), triangle.mrf$edges, triangle.mrf$node.par, triangle.mrf$edge.par, f)

# "Model matrix":
M <- t(sapply(1:nrow(config.mat), function(xx){phi.features(config = config.mat[xx,], triangle.mrf$edges, triangle.mrf$node.par, triangle.mrf$edge.par, f)}))
M
colSums(M)
#

# Compute sufficient statistics for all configurations:
s <- mrf.stat(triangle.mrf, instances = as.matrix(config.mat))
s

# Compute expexted number of features:
num.nodes <- length(gp@nodes)
E.phi     <- s/2^num.nodes
E.phi

set.seed(1)
N <- 5
samps <- sample.exact(triangle.mrf, size = N)
samps

s.emp <- mrf.stat(triangle.mrf, instances = samps)
s.emp

E.hat.phi <- s.emp/N
E.hat.phi
