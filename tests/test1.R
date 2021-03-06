# Initalize an mrf-object:
library(CRFutil)
library(Rgraphviz)


# Graph formula:
#grphf <- ~AA:BB:CC:DD
#grphf <- ~1:2 + 2:3 + 3:4 + 4:1
grphf <- ~A:B + B:C + C:D + D:A
#grphf <- ~A:B + C:B:A + D

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
Psi4 <- c(0.9,  0.1) *10

Psi12 <-
  rbind(
    c(30, 5),
    c(1, 10)
  )
Psi23 <-
  rbind(
    c(100, 1),
    c(1, 100)
  )
Psi34 <-
  rbind(
    c(1, 100),
    c(100, 1)
  )
Psi41 <-
  rbind(
    c(100, 1),
    c(1, 100)
  )

mc$node.pot[1,] <- Psi1
mc$node.pot[2,] <- Psi2
mc$node.pot[3,] <- Psi3
mc$node.pot[4,] <- Psi4

mc$edges # Check!
mc$edge.pot[[1]] <- Psi12
mc$edge.pot[[2]] <- Psi41
mc$edge.pot[[3]] <- Psi23
mc$edge.pot[[4]] <- Psi34

# Check again!
mc$node.pot
mc$edge.pot


#ar_new(names = "A", levels=list("A"=c("s1","s2")), values=c(1,2))
junk <- make.gRbase.potentials(mc, node.names = c("A","B","C","D"))
junk <- make.gRbase.potentials(mc, node.names = c("A","B","C","D"), state.nmes = c("+1","-1"))
junk <- make.gRbase.potentials(mc, node.names = gp@nodes, state.nmes = c("+1","-1"))
dimnames(junk$edge.potentials[[4]])
junk
log(junk$node.potentials[[1]])

# Energy func that takes one and two body energies
# Set and assign to MyEnergy??
s1<-1
s2<-2
st.sp <- expand.grid(c(s1,s2),c(s1,s2),c(s1,s2),c(s1,s2))
st.sp
pot.info <- make.gRbase.potentials(mc, node.names = gp@nodes)
pot.info
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}    # CRF states ??

energy.fun1(st.sp[1,], mc$edges, pot.info$node.energies, pot.info$edge.energies, f0)
E.vec2 <- sapply(1:nrow(st.sp), function(xx){energy.fun1(st.sp[xx,], mc$edges, pot.info$node.energies, pot.info$edge.energies, f0)})
E.vec2
exp(E.vec2) # prodPot's (P-tilde's)
sum(exp(E.vec2)) # Z

exp(E.vec2)/sum(exp(E.vec2)) # Prs

log(sum(exp(E.vec2))) # logZ us
infer.exact(mc)$logZ  # The correct logZ
infer.junction(mc)$logZ
infer.lbp(mc)$logZ
infer.trbp(mc)$logZ
infer.tree(mc)$logZ

exp(E.vec2)/exp(infer.exact(mc)$logZ)
exp(E.vec2)/exp(infer.junction(mc)$logZ)
exp(E.vec2)/exp(infer.lbp(mc)$logZ)
exp(E.vec2)/exp(infer.trbp(mc)$logZ)
exp(E.vec2)/exp(infer.tree(mc)$logZ)

sum(exp(E.vec2)/exp(infer.exact(mc)$logZ))
sum(exp(E.vec2)/exp(infer.junction(mc)$logZ))
sum(exp(E.vec2)/exp(infer.lbp(mc)$logZ))
sum(exp(E.vec2)/exp(infer.trbp(mc)$logZ))
sum(exp(E.vec2)/exp(infer.tree(mc)$logZ))

pot.info <- make.gRbase.potentials(mc, node.names = gp@nodes)
pot.info <- make.gRbase.potentials(mc, node.names = gp@nodes, state.nmes = c("+1","-1"))
pot.info

gR.dist.info <- distribution.from.potentials(pot.info$node.potentials, pot.info$edge.potentials)
as.data.frame(as.table(gR.dist.info$state.probs))

en.dist.info <- distribution.from.energies(st.sp, mc$edges, pot.info$node.energies, pot.info$edge.energies, energy, f0)
cbind(st.sp,en.dist.info$state.probs)


class(gR.dist.info$state.probs)
class(en.dist.info$state.probs)
