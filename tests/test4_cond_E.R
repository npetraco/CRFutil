# Initalize an mrf-object:
library(CRFutil)
library(Rgraphviz)


# Graph formula:
grphf <- ~A:B + A:C + B:C + C:D

# Check the graph:
gp <- ug(grphf, result = "graph")
dev.off()
plot(gp)
dev.off()
iplot(gp)

adj <- ug(grphf, result="matrix")
adj
n.states <- 2

mc <- make.crf(adj, n.states)
Psi1 <- c(0.25, 0.75)*4
Psi2 <- c(0.9,  0.1) *10
Psi3 <- c(0.25, 0.75)*4
Psi4 <- c(0.21, 0.87)*4

Psi12 <-
  rbind(
    c(30, 5),
    c(1, 10)
  )
Psi13 <-
  rbind(
    c(100, 1),
    c(1, 100)
  )
Psi23 <-
  rbind(
    c(1, 100),
    c(100, 1)
  )
Psi34 <-
  rbind(
    c(25, 10),
    c(6, 12)
  )


mc$node.pot[1,] <- Psi1
mc$node.pot[2,] <- Psi2
mc$node.pot[3,] <- Psi3
mc$node.pot[4,] <- Psi4

mc$edges # Check!
mc$edge.pot[[1]] <- Psi12
mc$edge.pot[[2]] <- Psi13
mc$edge.pot[[3]] <- Psi23
mc$edge.pot[[4]] <- Psi34

# Energy func that takes one and two body energies
# Set and assign to MyEnergy??
s1<-1
s2<-2
f0 <- function(y){ as.numeric(c((y==s1),(y==s2)))}

fake.sample <- rbind(
  c(1,1,1,2),
  c(1,2,1,1),
  c(2,2,2,2),
  c(1,2,1,2),
  c(2,1,1,2),
  c(2,1,2,2),
  c(1,1,1,1)
)


npe.info <- make.gRbase.potentials(mc, node.names = gp@nodes, state.nmes = c(1,2))
npe.info
mc$adj.nodes
mc$edges
conditional.config.energy(fake.sample[2,],
                          condition.element.number=3,
                          adj.node.list = mc$adj.nodes,
                          edge.mat = mc$edges,
                          one.lgp = npe.info$node.energies,
                          two.lgp = npe.info$edge.energies,
                          ff = f0)


samp.num <- 2
idx.interest <- 3
fake.sample[samp.num,]
e3 <- Eone(fake.sample[samp.num, idx.interest], npe.info$node.energies[[idx.interest]], f0)
e3
#31 = 13
fake.sample[samp.num, 1]
mc$edges
e13 <- Etwo(yA = fake.sample[samp.num, 1], yB = fake.sample[samp.num, idx.interest], wAB = npe.info$edge.energies[[2]], ff = f0)
e13
#32 = 23
fake.sample[samp.num, 2]
mc$edges
e23 <- Etwo(yA = fake.sample[samp.num, 2], yB = fake.sample[samp.num, idx.interest], wAB = npe.info$edge.energies[[3]], ff = f0)
e23
#34
fake.sample[samp.num, 4]
mc$edges
e34 <- Etwo(yA = fake.sample[samp.num, idx.interest], yB = fake.sample[samp.num, 4], wAB = npe.info$edge.energies[[4]], ff = f0)
e34

e3 + e13 + e23 + e34
exp(e3 + e13 + e23 + e34)/(1+exp(e3 + e13 + e23 + e34))

fake.sample[2,]
e.dn <- conditional.config.energy(
  c(1,2,2,1),
  condition.element.number=3,
  adj.node.list = mc$adj.nodes,
  edge.mat = mc$edges,
  one.lgp = npe.info$node.energies,
  two.lgp = npe.info$edge.energies,
  ff = f0)
e.up <- conditional.config.energy(
  c(1,2,1,1),
  condition.element.number=3,
  adj.node.list = mc$adj.nodes,
  edge.mat = mc$edges,
  one.lgp = npe.info$node.energies,
  two.lgp = npe.info$edge.energies,
  ff = f0)

e.up
e.dn

exp(e.up)/(1+exp(e.up))
exp(e.dn)/(1+exp(e.dn))
