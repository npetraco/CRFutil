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
num.samps <- 1000
#set.seed(1)
samps <- sample.exact(knm, num.samps)
mrf.sample.plot(samps)
#write.csv(samps, file = "/Users/npetraco/latex/papers/dust/steph_diss/CRFutil/tests/regression_tests/KNM_set1.csv")

knm.bel <- infer.exact(knm)
knm.bel$node.bel
knm.bel$edge.bel

# samps
# y1 <- samps[,1]
# y1[which(y==2)]<- 0
# glm(y ~ samps[,2] + samps[,3], family=binomial(link="logit"))

X <- samps
X1 <- X[,1]
X1[which(X1==2)]<- 0
X1
X2 <- X[,2]
X2[which(X2==2)]<- 0
X2
X3 <- X[,3]
X3[which(X3==2)]<- 0
X3

M1 <- glm(X1 ~ X2 + X3, family=binomial(link="logit"))
M2 <- glm(X2 ~ X1 + X3, family=binomial(link="logit"))
M3 <- glm(X3 ~ X1 + X2, family=binomial(link="logit"))

summary(M1)
coef(M1)[1] #     theta 1
coef(M1)[2] # 1-2 theta 4
coef(M1)[3] # 1-3 theta 5
# (Intercept)   -2.537 psl$node.par[1,1,] <- 1            theta 1
# X2             2.271 1-2 psl$edge.par[[1]][1,1,1] <- 4, theta 4
# X3             1.851 1-3 psl$edge.par[[2]][1,1,1] <- 5  theta 5
summary(M2)
coef(M2)[1] #     theta 2
coef(M2)[2] # 2-1 theta 4 *
coef(M2)[3] # 2-3 theta 6
# (Intercept)  -0.3877 psl$node.par[2,1,] <- 2            theta 2
# X1            2.2711 2-1 theta 4 too ??                 theta 4 *
# X3           -1.6694 2-3 psl$edge.par[[3]][1,1,1] <- 6  theta 6
summary(M3)
coef(M3)[1] #     theta 3
coef(M3)[2] # 3-1 theta 5 *
coef(M3)[3] # 3-2 theta 6 *
# (Intercept)  -1.1428 psl$node.par[3,1,] <- 3            theta 3
# X1            1.8510 3-1 theta 5??                      theta 5 *
# X2           -1.6694 3-2 theta 6??                      theta 6 *

#c(-2.537, -0.3877, -1.1428, 2.271, 1.851, -1.6694)
# knm$par
# -0.3510298  0.6039658 -1.2618153  0.6731633  0.4364151  0.1743194

#psl$par <- c(-2.537, -0.3877, -1.1428, 2.271, 1.851, -1.6694)
coef(M1)[1] #     theta 1
coef(M2)[1] #     theta 2
coef(M3)[1] #     theta 3
coef(M1)[2] # 1-2 theta 4
coef(M2)[2] # 2-1 theta 4 *
coef(M1)[3] # 1-3 theta 5
coef(M3)[2] # 3-1 theta 5 *
coef(M2)[3] # 2-3 theta 6
coef(M3)[3] # 3-2 theta 6 *

psl$par <- c(
  coef(M1)[1], #     theta 1
  coef(M2)[1], #     theta 2
  coef(M3)[1], #     theta 3
  coef(M1)[2], # 1-2 theta 4
  coef(M1)[3], # 1-3 theta 5
  coef(M2)[3]) # 2-3 theta 6


out.pot2 <- make.pots(parms = psl$par,  crf = psl,  rescaleQ = T, replaceQ = T)
psl$node.pot
psl$edge.pot

psl.bel <- infer.exact(psl)
psl.bel$node.bel
knm.bel$node.bel

psl.bel$edge.bel[[1]]
knm.bel$edge.bel[[1]]
psl.bel$edge.bel[[2]]
knm.bel$edge.bel[[2]]
psl.bel$edge.bel[[3]]
knm.bel$edge.bel[[3]]
