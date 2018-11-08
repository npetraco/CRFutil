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

configs <- rbind(
  c(1,1,1), #X1
  c(2,1,1), #X2
  c(1,2,1), #X3
  c(2,2,1), #X4
  c(1,1,2), #X5
  c(2,1,2), #X6
  c(1,2,2), #X7
  c(2,2,2)  #X8
)

X  <- configs[8,]
Xc <- complement.at.idx(configuration = X, complement.index = 3)
symbolic.conditional.energy(config = X,  condition.element.number = 3, crf = psl, ff = f0, printQ = F)
symbolic.conditional.energy(config = Xc, condition.element.number = 3, crf = psl, ff = f0, printQ = F)

#phi.component(config = X, i=3, node.par=psl$node.par, ff=f0)
#phi.component(config = X, i=1, j=3, edge.par = psl$edge.par, edge.mat = psl$edges, ff = f0)
ec1  <- symbolic.conditional.energy(config = X,  condition.element.number = 1, crf = psl, ff = f0, printQ = F, format = "conditional.phi")
ecc1 <- symbolic.conditional.energy(config = complement.at.idx(configuration = X, complement.index = 1),  condition.element.number = 1, crf = psl, ff = f0, printQ = F, format = "conditional.phi")
ec1-ecc1

ec2  <- symbolic.conditional.energy(config = X,  condition.element.number = 2, crf = psl, ff = f0, printQ = F, format = "conditional.phi")
ecc2 <- symbolic.conditional.energy(config = complement.at.idx(configuration = X, complement.index = 2),  condition.element.number = 2, crf = psl, ff = f0, printQ = F, format = "conditional.phi")
ec2-ecc2

ec3  <- symbolic.conditional.energy(config = X,  condition.element.number = 3, crf = psl, ff = f0, printQ = F, format = "conditional.phi")
ecc3 <- symbolic.conditional.energy(config = complement.at.idx(configuration = X, complement.index = 3),  condition.element.number = 3, crf = psl, ff = f0, printQ = F, format = "conditional.phi")
ec3-ecc3

(ec1-ecc1) + (ec2-ecc2) + (ec3-ecc3)

# Delta-alpha matrix:
MX <- array(NA,c(nrow(configs)*psl$n.nodes, psl$n.par))

count <- 1
for(i in 1:nrow(configs)) {
  for(j in 1:psl$n.nodes) {
    ec  <- symbolic.conditional.energy(config = configs[i,],  condition.element.number = j, crf = psl, ff = f0, printQ = F, format = "conditional.phi")
    ecc <- symbolic.conditional.energy(config = complement.at.idx(configuration = configs[i,], complement.index = j),  condition.element.number = j, crf = psl, ff = f0, printQ = F, format = "conditional.phi")
    print(paste("sample X:", i, "Node:", j) )
    #print(ec)
    #print(ecc)
    print(ec-ecc)
    MX[count,] <- ec-ecc
    count <- count + 1
  }
}
MX

nrow(configs)*psl$n.nodes
