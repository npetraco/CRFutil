library(CRFutil)

grphf <- ~A:B + B:C + C:D + D:A
adj   <- ug(grphf, result="matrix")

n.states <- 2
com <- make.crf(adj, n.states)
com <- make.features(com)
com <- make.par(com, 12)

# Fill the parameter index matrices:
for(i in 1:nrow(com$node.par)){
  com$node.par[i,1,1] <- i
}
com$edge.par[[1]][1,1,1] <- 5
com$edge.par[[1]][2,2,1] <- 6
com$edge.par[[2]][1,1,1] <- 7
com$edge.par[[2]][2,2,1] <- 8
com$edge.par[[3]][1,1,1] <- 9
com$edge.par[[3]][2,2,1] <- 10
com$edge.par[[4]][1,1,1] <- 11
com$edge.par[[4]][2,2,1] <- 12

com$node.par
com$edge.par
com$edges
nodes2params.list(com, storeQ = F)


# Fill the parameter index matrices:

com$node.par[1,1,1] <- 1
com$node.par[2,1,1] <- 2
com$node.par[3,1,1] <- 3
com$node.par[4,1,1] <- 4

com$edge.par[[1]][1,1,1] <- 5
com$edge.par[[1]][2,2,1] <- 6
com$edge.par[[1]][2,1,1] <- 7

com$edge.par[[2]][1,1,1] <- 8
com$edge.par[[2]][2,2,1] <- 9
com$edge.par[[2]][2,1,1] <- 10

com$edge.par[[3]][1,1,1] <- 11
com$edge.par[[3]][2,2,1] <- 12
com$edge.par[[3]][2,1,1] <- 13

com$edge.par[[4]][1,1,1] <- 14
com$edge.par[[4]][2,2,1] <- 15
com$edge.par[[4]][2,1,1] <- 16

com$node.par
com$edge.par
com$edges
nodes2params.list(com, storeQ = F)


#---------------------------------------
grphf <- ~A:B + B:C + C:D + D:A + A:C
adj   <- ug(grphf, result="matrix")
iplot(ug(grphf))
adj

n.states <- 2
com <- make.crf(adj, n.states)
com <- make.features(com)
com <- make.par(com, 23)

com$node.par[1,,1] <- c(1,2)
com$node.par[2,,1] <- c(3,4)
com$node.par[3,,1] <- c(5,6)
com$node.par[4,,1] <- c(7,8)
com$node.par

com$edge.par[[1]][1,1,1] <- 9
com$edge.par[[1]][2,2,1] <- 10
com$edge.par[[1]][2,1,1] <- 11

com$edge.par[[2]][1,1,1] <- 12
com$edge.par[[2]][2,2,1] <- 13
com$edge.par[[2]][2,1,1] <- 14

com$edge.par[[3]][1,1,1] <- 15
com$edge.par[[3]][2,2,1] <- 16
com$edge.par[[3]][2,1,1] <- 17

com$edge.par[[4]][1,1,1] <- 18
com$edge.par[[4]][2,2,1] <- 19
com$edge.par[[4]][2,1,1] <- 20

com$edge.par[[5]][1,1,1] <- 21
com$edge.par[[5]][2,2,1] <- 22
com$edge.par[[5]][2,1,1] <- 23


com$node.par
com$edge.par
com$edges
nodes2params.list(com, storeQ = T)
params2nodes.list(com, storeQ = T)
com$pars2nodes
com$nodes2pars

