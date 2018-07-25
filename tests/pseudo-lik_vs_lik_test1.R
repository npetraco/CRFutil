library(CRFutil)

grphf <- ~A:B:C
adj   <- ug(grphf, result="matrix")
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}

n.states <- 2
known.model <- make.crf(adj, n.states)

# True node pots:
PsiA <- c(2,1)
PsiB <- c(1,3)
PsiC <- c(3,1)

# True edge pots:
PsiAB <-
  rbind(
    c(3, 6.1),
    c(6.1, 3.6)
  )

PsiBC <-
  rbind(
    c(2.5, 3.1),
    c(3.1, 2)
  )

PsiAC <-
  rbind(
    c(4, 1.1),
    c(1.1, 4.3)
  )

known.model$node.pot[1,]  <- PsiA
known.model$node.pot[2,]  <- PsiB
known.model$node.pot[2,]  <- PsiC

known.model$edge.pot[[1]] <- PsiAB
known.model$edge.pot[[2]] <- PsiAC
known.model$edge.pot[[3]] <- PsiBC

# So now sample from the model as if we obtained an experimental sample:
num.samps <- 25
set.seed(1)
samps <- sample.exact(known.model, num.samps)
mrf.sample.plot(samps)


library(rgl)

# Instantiate a two node two parameter model to fit:
fit <- make.crf(adj, n.states)
fit <- make.features(fit)
fit <- make.par(fit, 2)
fit$node.par[1,1,] <- 1
fit$node.par[2,1,] <- 2
fit$node.par[3,1,] <- 3
fit$node.par

fit$edge.par[[1]][1,1,1] <- 4
fit$edge.par[[1]][2,2,1] <- 5
fit$edge.par[[2]][1,1,1] <- 6
fit$edge.par[[2]][2,2,1] <- 7
fit$edge.par[[3]][1,1,1] <- 8
fit$edge.par[[3]][2,2,1] <- 9
fit$edge.par
n2p <- nodes2params.list2(fit, storeQ = T)
n2p
fit$adj.nodes


X.all <- expand.grid(c(1,2),c(1,2),c(1,2))
X.all
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}

# Log Pseudo-Likelihood:
t1 <- fit$node.par[1,,1]
t2 <- fit$node.par[2,,1]
t3 <- fit$node.par[3,,1]
w12 <- fit$edge.par[[1]][,,]
w13 <- fit$edge.par[[2]][,,]
w23 <- fit$edge.par[[3]][,,]

X <- c(1,1,1)
# X1
c( f0(X[1])%*%t1, f0(X[1])%*%w12%*%f0(X[2]), f0(X[1])%*%w13%*%f0(X[3]) )
# X2
c( f0(X[2])%*%t2, f0(X[2])%*%w12%*%f0(X[1]), f0(X[2])%*%w23%*%f0(X[3]) )
# X3
c( f0(X[3])%*%t3, f0(X[3])%*%w23%*%f0(X[2]), f0(X[3])%*%w13%*%f0(X[1]) )

# Complements:
X <- X.all[8,]
X
c(
  row.match(complement.at.idx(X, complement.index = 1), table = X.all),
  row.match(complement.at.idx(X, complement.index = 2), table = X.all),
  row.match(complement.at.idx(X, complement.index = 3), table = X.all) )

complement.at.idx(X, complement.index = 1)
complement.at.idx(X, complement.index = 2)
complement.at.idx(X, complement.index = 3)


# Log Likelihood:
X <- c(1,1,1)
phi.X <- phi.features(X, edges.mat = fit$edges, node.par = fit$node.par, edge.par = fit$edge.par, ff = f0)
theta <- 1:9
theta * phi.X
