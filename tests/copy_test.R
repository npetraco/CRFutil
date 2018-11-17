library(CRFutil)

# Fully connected
grphf <- ~1:2+1:3+2:3
adj <- ug(grphf, result="matrix")
# Check the graph:
gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)

f0 <- function(y){ as.numeric(c((y==1),(y==2)))}
n.states <- 2


knm <- make.crf(adj, n.states)
knm <- make.features(knm)
knm <- make.par(knm, 6)
knm$node.par[1,1,] <- 1
knm$node.par[2,1,] <- 2
knm$node.par[3,1,] <- 3
knm$edge.par[[1]][1,1,1] <- 4
knm$edge.par[[1]][2,2,1] <- 4
knm$edge.par[[2]][1,1,1] <- 5
knm$edge.par[[2]][2,2,1] <- 5
knm$edge.par[[3]][1,1,1] <- 6
knm$edge.par[[3]][2,2,1] <- 6

set.seed(6)
knm$par <- runif(6,-1.5,1.1)
out.pot <- make.pots(parms = knm$par,  crf = knm,  rescaleQ = T, replaceQ = T)

#----- Deep copy an existing model:
a.cpy.md <- copy.crf(knm, plotQ = T)
a.cpy.md$adj.edges
knm$adj.edges

a.cpy.md$adj.nodes
knm$adj.nodes

a.cpy.md$edge.par
knm$edge.par

a.cpy.md$edge.pot
knm$edge.pot

a.cpy.md$edges
knm$edges

a.cpy.md$gradient
knm$gradient

a.cpy.md$max.state
knm$max.state

a.cpy.md$n.adj
knm$n.adj

a.cpy.md$n.edges
knm$n.edges

a.cpy.md$n.ef
knm$n.ef

a.cpy.md$n.nf
knm$n.nf

a.cpy.md$n.nodes
knm$n.nodes

a.cpy.md$n.par
knm$n.par

a.cpy.md$n.states
knm$n.states

a.cpy.md$nll
knm$nll

a.cpy.md$node.par
knm$node.par

a.cpy.md$node.pot
knm$node.pot

a.cpy.md$par
knm$par

dump.crf(knm)

#----------- Instantiate an empty model:
# Adamantane field:
grphf <- ~1:2 + 2:3 + 3:4 + 4:5 + 5:6 + 6:1 + 1:7 + 7:8 + 8:9 + 9:6 + 10:8 + 10:3
adj <- ug(grphf, result="matrix")
# Check the graph:
gp <- ug(grphf, result = "graph")
dev.off()
iplot(gp)

grphf
adj

empty.modl <- make.empty.field(
  graph.eq             = NULL,
  adj.mat              = adj,
  parameterization.typ = "ising2",
  node.par             = NULL,
  edge.par             = NULL,
  plotQ                = FALSE)


dump.crf(empty.modl)
