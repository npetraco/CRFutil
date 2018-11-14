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

# Instantiate an empty model to fit:
psl <- make.crf(adj, n.states)
psl <- make.features(psl)
psl <- make.par(psl, 6)
psl$node.par[1,1,] <- 1
psl$node.par[2,1,] <- 2
psl$node.par[3,1,] <- 3

psl$edge.par[[1]][1,1,1] <- 4
psl$edge.par[[1]][2,2,1] <- 4
psl$edge.par[[2]][1,1,1] <- 5
psl$edge.par[[2]][2,2,1] <- 5
psl$edge.par[[3]][1,1,1] <- 6
psl$edge.par[[3]][2,2,1] <- 6
psl$edges


# All possible configs:
samps <- rbind(
  c(1,1,1), #X1
  c(2,1,1), #X2
  c(1,2,1), #X3
  c(2,2,1), #X4
  c(1,1,2), #X5
  c(2,1,2), #X6
  c(1,2,2), #X7
  c(2,2,2)  #X8
)

# Delta-alpha matrix:
MX      <- array(NA,c(nrow(samps)*psl$n.nodes, psl$n.par))
ec.mat  <- array(NA,c(nrow(samps)*psl$n.nodes, psl$n.par))
ecc.mat <- array(NA,c(nrow(samps)*psl$n.nodes, psl$n.par))


count <- 1
#for(i in 1:nrow(configs)) {
#  for(j in 1:psl$n.nodes) {
for(j in 1:psl$n.nodes) {
  for(i in 1:nrow(samps)) {

    # Convert to 1/0 node states
    X <- samps[i,]
    #X[which(X == 2)] <- 0

    # Make Xj = 1, Xjc = 0
    Xc    <- X
    X[j]  <- 1
    Xc[j] <- 2

    ec  <- symbolic.conditional.energy(config = X,  condition.element.number = j, crf = psl, ff = f0, printQ = F, format = "conditional.phi")
    ecc <- symbolic.conditional.energy(config = Xc, condition.element.number = j, crf = psl, ff = f0, printQ = F, format = "conditional.phi")

    txec  <- symbolic.conditional.energy(config = X,  condition.element.number = j, crf = psl, ff = f0, printQ = F, format = "tex")
    txecc <- symbolic.conditional.energy(config = Xc, condition.element.number = j, crf = psl, ff = f0, printQ = F, format = "tex")

    print(paste("sample X:", i, "Node:", j) )
    print(X)
    print(Xc)
    print(ec)
    print(ecc)
    print(txec)
    print(txecc)
    print(ec-ecc)
    print("=====================")
    ec.mat[count,]  <- ec
    ecc.mat[count,] <- ecc
    MX[count,]      <- ec-ecc
    count           <- count + 1
  }
}
dim(MX)
ec.mat
ecc.mat
MX
