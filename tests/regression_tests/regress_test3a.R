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

set.seed(6)
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
    print(paste("sample X:", i, "Node:", j) )
    print(X)
    print(Xc)
    print(ec)
    print(ecc)
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

cbind(MX[1:25, c(1,4,5)], samps)
cbind(MX[26:50,c(2,4,6)], samps)
cbind(MX[51:75,c(3,5,6)], samps)

y <-c(samps[,1], samps[,2], samps[,3])
y
y[which(y==2)] <- 0
y

Ma <- glm(y ~ MX[,1] + MX[,2] + MX[,3] + MX[,4] + MX[,5] + MX[,6] - 1, family=binomial(link="logit"))
summary(Ma)

MX1 <- MX[1:25, c(1,4,5)]
y1  <- y[1:25]
M1 <- glm(y1 ~ MX1[,1] + MX1[,2] + MX1[,3] - 1, family=binomial(link="logit"))
M1a <- glm(y1 ~ MX1[,2] + MX1[,3], family=binomial(link="logit"))
summary(M1)
summary(M1a)

coef(M1)
coef(M1a)
coef(Ma)[c(1,4,5)]



MX2 <- MX[26:50,c(2,4,6)]
y2  <- y[26:50]
M2 <- glm(y2 ~ MX2[,1] + MX2[,2] + MX2[,3] - 1, family=binomial(link="logit"))
M2a <- glm(y2 ~ MX2[,2] + MX2[,3], family=binomial(link="logit"))
summary(M2)
summary(M2a)

coef(M2)
coef(M2a)
coef(Ma)[c(2,4,6)]



MX3 <- MX[51:75,c(3,5,6)]
y3  <- y[51:75]
M3 <- glm(y3 ~ MX3[,1] + MX3[,2] + MX3[,3] - 1, family=binomial(link="logit"))
M3a <- glm(y3 ~ MX3[,2] + MX3[,3], family=binomial(link="logit"))
summary(M3)
summary(M3a)

coef(M3)
coef(M3a)
coef(Ma)[c(3,5,6)]



coef(Ma)
knm$par

psl$par <- as.numeric(coef(Ma))

out.pot2 <- make.pots(parms = psl$par,  crf = psl,  rescaleQ = T, replaceQ = T)
psl$node.pot
psl$edge.pot

# Node and edge beliefs:
psl.bel <- infer.exact(psl)
psl.bel$node.bel
knm.bel$node.bel

psl.bel$edge.bel[[1]]
knm.bel$edge.bel[[1]]
psl.bel$edge.bel[[2]]
knm.bel$edge.bel[[2]]
psl.bel$edge.bel[[3]]
knm.bel$edge.bel[[3]]

# Configuration probabilities:
pot.info <- make.gRbase.potentials(psl, node.names = gp@nodes, state.nmes = c("1","2"))
pot.info

gR.dist.info    <- distribution.from.potentials(pot.info$node.potentials, pot.info$edge.potentials)
logZ            <- gR.dist.info$logZ
joint.dist.info <- as.data.frame(as.table(gR.dist.info$state.probs))
joint.dist.info

pot.info.knm <- make.gRbase.potentials(knm, node.names = gp@nodes, state.nmes = c("1","2"))
pot.info.knm

gR.dist.info.knm    <- distribution.from.potentials(pot.info.knm$node.potentials, pot.info.knm$edge.potentials)
logZ.knm            <- gR.dist.info.knm$logZ
joint.dist.info.knm <- as.data.frame(as.table(gR.dist.info.knm$state.probs))
joint.dist.info.knm

cbind(joint.dist.info, joint.dist.info.knm[,4])

