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

knm$par <- runif(6,-1.5,1.1)
knm$par # "true" theta
#knm$par <-
out.pot <- make.pots(parms = knm$par,  crf = knm,  rescaleQ = T, replaceQ = T)
knm$node.pot
knm$edge.pot

# So now sample from the model as if we obtained an experimental sample:
num.samps <- 25
#set.seed(1)
samps <- sample.exact(knm, num.samps)
mrf.sample.plot(samps)
#write.csv(samps, file = "/Users/npetraco/latex/papers/dust/steph_diss/CRFutil/tests/regression_tests/KNM_set1.csv")

knm.bel <- infer.exact(knm)
knm.bel$node.bel
knm.bel$edge.bel


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

# Delta-alpha matrix:
MX <- array(NA,c(nrow(samps)*psl$n.nodes, psl$n.par))

count <- 1
#for(i in 1:nrow(configs)) {
#  for(j in 1:psl$n.nodes) {
for(j in 1:psl$n.nodes) {
  for(i in 1:nrow(samps)) {
    ec  <- symbolic.conditional.energy(config = samps[i,],  condition.element.number = j, crf = psl, ff = f0, printQ = F, format = "conditional.phi")
    ecc <- symbolic.conditional.energy(config = complement.at.idx(configuration = samps[i,], complement.index = j),  condition.element.number = j, crf = psl, ff = f0, printQ = F, format = "conditional.phi")
    print(paste("sample X:", i, "Node:", j) )
    #print(ec)
    #print(ecc)
    print(ec-ecc)
    MX[count,] <- ec-ecc
    count <- count + 1
  }
}
dim(MX)

y <-c(samps[,1], samps[,2], samps[,3])
y
y[which(y==2)] <- 0
y

M1 <- glm(y ~ MX[,1] + MX[,2] + MX[,3] + MX[,4] + MX[,5] + MX[,6] - 1, family=binomial(link="logit"))
summary(M1)
