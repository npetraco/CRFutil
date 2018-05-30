# Initalize an mrf-object:
library(CRFutil)
library(Rgraphviz)


# Graph formula:
grphf <- ~A:B + B:C + C:D

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
# These are what CRF takes as inputs/fits. The "potentials".
Psi1 <- c(0.25, 0.75)*4
Psi2 <- c(0.9,  0.1) *10
Psi3 <- c(0.25, 0.75)*4
Psi4 <- c(0.9,  0.1) *10

Psi12 <-
  6*rbind(c(2/6, 1/6),
          c(1/6, 2/6))
Psi23 <-
  6*rbind(c(2/6, 1/6),
          c(1/6, 2/6))
Psi34 <-
  6*rbind(c(2/6, 1/6),
          c(1/6, 2/6))


mc$node.pot[1,] <- Psi1
mc$node.pot[2,] <- Psi2
mc$node.pot[3,] <- Psi3
mc$node.pot[4,] <- Psi4

mc$edges # Check!
mc$edge.pot[[1]] <- Psi12
mc$edge.pot[[2]] <- Psi23
mc$edge.pot[[3]] <- Psi34

# Check again!
mc$node.pot
mc$edge.pot

s1<-1
s2<-2
st.sp <- expand.grid(c(s1,s2),c(s1,s2),c(s1,s2),c(s1,s2))
st.sp
pot.info <- make.gRbase.potentials(mc, node.names = gp@nodes)
f0 <- function(y){ as.numeric(c((y==1),(y==2)))}    # Feature function

#
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

