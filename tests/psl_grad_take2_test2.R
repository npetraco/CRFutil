library(CRFutil)

# Graph formula for Slayer field:
grphf <- ~1:2+1:3+1:4

# Check the graph:
gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)

# Adjacenty matrix:
adj <- ug(grphf, result="matrix")

# Define states and feature function:
s1 <- 1
s2 <- 2
f0 <- function(y){ as.numeric(c((y==s1),(y==s2))) }

n.states <- 2
fitc <- make.crf(adj, n.states)
fitc <- make.features(fitc)
fitc <- make.par(fitc, 6)

fitc$edges

# Parameterization:
fitc$node.par
fitc$node.par[1,1,] <- 1 #1
fitc$node.par[2,1,] <- 2 #2
fitc$node.par[3,1,] <- 1 #1
fitc$node.par[4,1,] <- 3 #3

fitc$edge.par
fitc$edge.par[[1]][1,1,1] <- 1 #1
fitc$edge.par[[1]][2,2,1] <- 4 #4

fitc$edge.par[[2]][1,1,1] <- 1 #1
fitc$edge.par[[2]][2,2,1] <- 5 #5

fitc$edge.par[[3]][1,1,1] <- 6 #6
fitc$edge.par[[3]][2,2,1] <- 1 #1

fitc$par

X <- rbind(
  c(1,1,1,1),
  c(1,2,1,1),
  c(1,1,2,1),
  c(1,1,1,2),
  c(2,1,1,1),
  c(2,2,1,1),
  c(2,1,2,1),
  c(2,1,1,2),
  c(1,2,2,1),
  c(1,2,1,2),
  c(1,1,2,2),
  c(1,2,2,2),
  c(2,1,2,2),
  c(2,2,1,2),
  c(2,2,2,1),
  c(2,2,2,2)
)
X

conditional.energy.gradient(config = X[1,], condition.element.number = 4, crf = fitc, ff = f0, printQ=T)

# E(X_i|X/X_i)  i = 1, X=(1,1,1,1)
#ell contribution from node i=1 \theta_{\ell[i]}
get.par.idx(config =X[1,], i=1, node.par=fitc$node.par, ff=f0)

#m cotribution from edge i=1,j=2 \theta_{m[i~j]}
fitc$edges
get.par.idx(config   = X[1,],
            i        = 1,
            j        = 2,
            edge.par = fitc$edge.par,
            edge.mat = fitc$edges,
            ff       = f0)

#m cotribution from edge i=1,j=3 \theta_{m[i~j]}
fitc$edges
get.par.idx(config   = X[1,],
            i        = 1,
            j        = 3,
            edge.par = fitc$edge.par,
            edge.mat = fitc$edges,
            ff       = f0)

#m cotribution from edge i=1,j=4 \theta_{m[i~j]}
fitc$edges
get.par.idx(config   = X[1,],
            i        = 1,
            j        = 4,
            edge.par = fitc$edge.par,
            edge.mat = fitc$edges,
            ff       = f0)

# So E(X_1|X/X_1) = theta_1 + theta_1 + theta_1 + theta_6
#
X <- c(1,1,1,1)
conditional.energy.gradient(config = X, condition.element.number = 1, crf = fitc, ff = f0)$conditional.grad
conditional.energy.gradient(config = X, condition.element.number = 2, crf = fitc, ff = f0)$conditional.grad
conditional.energy.gradient(config = X, condition.element.number = 3, crf = fitc, ff = f0)$conditional.grad
conditional.energy.gradient(config = X, condition.element.number = 4, crf = fitc, ff = f0)$conditional.grad

symbolic.conditional.energy(config = X, condition.element.number = 4, crf = fitc, ff = f0)
