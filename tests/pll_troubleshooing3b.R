# set.seed(1)
# tmp <- runif(2000, min = -10, max = 10)
# thetas <- cbind(tmp[1:1000], tmp[1001:2000])
# write.csv(x = thetas, file = paste0(getwd(), "/tests/rnd_thetas.csv"))

# Node numbers
nods <- c(1,2)
# Possible state configurations:
configs <- rbind(
  c(1,1),
  c(1,2),
  c(2,1),
  c(2,2)
)


en.info     <- array(NA,c(nrow(thetas)*nrow(configs)*length(nods),4))
en.info2    <- array(NA,c(nrow(thetas)*nrow(configs)*length(nods),4))
count <- 1
for(i in 1:nrow(thetas)) {
  for(j in 1:nrow(configs)) {
    for(k in 1:length(nods)) {

      theta    <- thetas[i,]
      X       <- configs[j,]
      node.num <- nods[k]

      # Original recipie:
      theta.reformated <- par2logpots(theta, fitc)
      en <- conditional.config.energy(
        config                   = X,
        condition.element.number = node.num,
        adj.node.list            = fitc$adj.nodes,
        edge.mat                 = fitc$edges,
        one.lgp                  = theta.reformated[[1]],
        two.lgp                  = theta.reformated[[2]],
        ff                       = f0,
        printQ                   = F)
      en.info[count,] <- c(i,j,k,en)

      # New recipie:
      #Compute phi for X:
      phi.X <- phi.features(
        config    = X,
        edges.mat = fitc$edges,
        node.par  = fitc$node.par,
        edge.par  = fitc$edge.par,
        ff        = f0
      )

      # Grab the parameters associated with conditional node (node.num)
      n2p <- nodes2params.list(fitc, storeQ = T)
      node.pars <- n2p[[node.num]]

      # Compute the energy: E(X_i|{\bf X}_i\slash X_i) using the dot product formulation
      fitc$par <- theta
      en2 <- fitc$par[node.pars] %*% phi.X[node.pars]
      en.info2[count,] <- c(i,j,k,en2)

      count <- count + 1

    }
  }
}

en.info[,4]
en.info2[,4]
en.info[,4] - en.info2[,4]

hist(en.info[,4])


#thetas = Import["C:\\Users\\Victor Lin\\Downloads\\rnd_thetas.csv"];
#config = {{1, 1}, {1, 2}, {2, 1}, {2, 2}};
tab <- c("theta1", "theta2", "(1,1)e1", "(1,1)e2", "(1,2)e1", "(1,2)e2", "(2,1)e1", "(2,1)e2", "(2,2)e1", "(2,2)e2")
r.all  <- NULL
r.all2 <- NULL
for(tnum in 1:nrow(thetas)){
  r  <- thetas[tnum,]
  r2 <- thetas[tnum,]
  for(nnum in 1:4) {
    for(num in 1:2) {

      theta    <- thetas[tnum,]
      X       <- configs[nnum,]
      node.num <- nods[num]

      # Original recipie:
      theta.reformated <- par2logpots(theta, fitc)
      en <- conditional.config.energy(
        config                   = X,
        condition.element.number = node.num,
        adj.node.list            = fitc$adj.nodes,
        edge.mat                 = fitc$edges,
        one.lgp                  = theta.reformated[[1]],
        two.lgp                  = theta.reformated[[2]],
        ff                       = f0,
        printQ                   = F)
      #en.info[count,] <- c(i,j,k,en)
      #count <- count + 1

      # New recipie:
      #Compute phi for X:
      phi.X <- phi.features(
        config    = X,
        edges.mat = fitc$edges,
        node.par  = fitc$node.par,
        edge.par  = fitc$edge.par,
        ff        = f0
      )

      # Grab the parameters associated with conditional node (node.num)
      n2p <- nodes2params.list(fitc, storeQ = T)
      node.pars <- n2p[[node.num]]

      # Compute the energy: E(X_i|{\bf X}_i\slash X_i) using the dot product formulation
      fitc$par <- theta
      en2 <- fitc$par[node.pars] %*% phi.X[node.pars]

      r  <- c(r,  en)
      r2 <- c(r2, en2)
    }
  }
  r.all  <- rbind(r.all, r)
  r.all2 <- rbind(r.all2,r2)
}
colnames(r.all)  <- tab
colnames(r.all2) <- tab
r.all
r.all2

r.all.vic <- read.csv(paste0(getwd(), "/tests/rnd_thetas_vic_output.csv"), header = T)
r.all[,3:10] - r.all.vic[,3:10]
max(r.all[,3:10] - r.all.vic[,3:10])
min(r.all[,3:10] - r.all.vic[,3:10])

r.all[,3:10] - r.all2[,3:10]
