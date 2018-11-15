library(CRFutil)


# Fully connected
grphf <- ~1:2 + 1:3 + 1:4 + 1:5 + 1:6 + 1:7 + 1:8 + 1:9 + 1:10 + 2:3 + 2:4 + 2:5 + 2:6 + 2:7 + 2:8 + 2:9 + 2:10 + 3:4 + 3:5 + 3:6 + 3:7 + 3:8 + 3:9 + 3:10 + 4:5 + 4:6 + 4:7 + 4:8 + 4:9 + 4:10 + 5:6 + 5:7 + 5:8 + 5:9 + 5:10 + 6:7 + 6:8 + 6:9 + 6:10 + 7:8 + 7:9 + 7:10 + 8:9 + 8:10 + 9:10
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
knm <- make.par(knm, 55)
knm$n.edges

knm$node.par[1,1,] <- 1
knm$node.par[2,1,] <- 2
knm$node.par[3,1,] <- 3
knm$node.par[4,1,] <- 4
knm$node.par[5,1,] <- 5
knm$node.par[6,1,] <- 6
knm$node.par[7,1,] <- 7
knm$node.par[8,1,] <- 8
knm$node.par[9,1,] <- 9
knm$node.par[10,1,] <- 10
knm$edge.par[[1]][1,1,1] <- 11
knm$edge.par[[1]][2,2,1] <- 11
knm$edge.par[[2]][1,1,1] <- 12
knm$edge.par[[2]][2,2,1] <- 12
knm$edge.par[[3]][1,1,1] <- 13
knm$edge.par[[3]][2,2,1] <- 13
knm$edge.par[[4]][1,1,1] <- 14
knm$edge.par[[4]][2,2,1] <- 14
knm$edge.par[[5]][1,1,1] <- 15
knm$edge.par[[5]][2,2,1] <- 15
knm$edge.par[[6]][1,1,1] <- 16
knm$edge.par[[6]][2,2,1] <- 16
knm$edge.par[[7]][1,1,1] <- 17
knm$edge.par[[7]][2,2,1] <- 17
knm$edge.par[[8]][1,1,1] <- 18
knm$edge.par[[8]][2,2,1] <- 18
knm$edge.par[[9]][1,1,1] <- 19
knm$edge.par[[9]][2,2,1] <- 19
knm$edge.par[[10]][1,1,1] <- 20
knm$edge.par[[10]][2,2,1] <- 20

knm$edge.par[[11]][1,1,1] <- 21
knm$edge.par[[11]][2,2,1] <- 21
knm$edge.par[[12]][1,1,1] <- 22
knm$edge.par[[12]][2,2,1] <- 22
knm$edge.par[[13]][1,1,1] <- 23
knm$edge.par[[13]][2,2,1] <- 23
knm$edge.par[[14]][1,1,1] <- 24
knm$edge.par[[14]][2,2,1] <- 24
knm$edge.par[[15]][1,1,1] <- 25
knm$edge.par[[15]][2,2,1] <- 25
knm$edge.par[[16]][1,1,1] <- 26
knm$edge.par[[16]][2,2,1] <- 26
knm$edge.par[[17]][1,1,1] <- 27
knm$edge.par[[17]][2,2,1] <- 27
knm$edge.par[[18]][1,1,1] <- 28
knm$edge.par[[18]][2,2,1] <- 28
knm$edge.par[[19]][1,1,1] <- 29
knm$edge.par[[19]][2,2,1] <- 29
knm$edge.par[[20]][1,1,1] <- 30
knm$edge.par[[20]][2,2,1] <- 30

knm$edge.par[[21]][1,1,1] <- 31
knm$edge.par[[21]][2,2,1] <- 31
knm$edge.par[[22]][1,1,1] <- 32
knm$edge.par[[22]][2,2,1] <- 32
knm$edge.par[[23]][1,1,1] <- 33
knm$edge.par[[23]][2,2,1] <- 33
knm$edge.par[[24]][1,1,1] <- 34
knm$edge.par[[24]][2,2,1] <- 34
knm$edge.par[[25]][1,1,1] <- 35
knm$edge.par[[25]][2,2,1] <- 35
knm$edge.par[[26]][1,1,1] <- 36
knm$edge.par[[26]][2,2,1] <- 36
knm$edge.par[[27]][1,1,1] <- 37
knm$edge.par[[27]][2,2,1] <- 37
knm$edge.par[[28]][1,1,1] <- 38
knm$edge.par[[28]][2,2,1] <- 38
knm$edge.par[[29]][1,1,1] <- 39
knm$edge.par[[29]][2,2,1] <- 39
knm$edge.par[[30]][1,1,1] <- 40
knm$edge.par[[30]][2,2,1] <- 40

knm$edge.par[[31]][1,1,1] <- 41
knm$edge.par[[31]][2,2,1] <- 41
knm$edge.par[[32]][1,1,1] <- 42
knm$edge.par[[32]][2,2,1] <- 42
knm$edge.par[[33]][1,1,1] <- 43
knm$edge.par[[33]][2,2,1] <- 43
knm$edge.par[[34]][1,1,1] <- 44
knm$edge.par[[34]][2,2,1] <- 44
knm$edge.par[[35]][1,1,1] <- 45
knm$edge.par[[35]][2,2,1] <- 45
knm$edge.par[[36]][1,1,1] <- 46
knm$edge.par[[36]][2,2,1] <- 46
knm$edge.par[[37]][1,1,1] <- 47
knm$edge.par[[37]][2,2,1] <- 47
knm$edge.par[[38]][1,1,1] <- 48
knm$edge.par[[38]][2,2,1] <- 48
knm$edge.par[[39]][1,1,1] <- 49
knm$edge.par[[39]][2,2,1] <- 49
knm$edge.par[[40]][1,1,1] <- 50
knm$edge.par[[40]][2,2,1] <- 50

knm$edge.par[[41]][1,1,1] <- 51
knm$edge.par[[41]][2,2,1] <- 51
knm$edge.par[[42]][1,1,1] <- 52
knm$edge.par[[42]][2,2,1] <- 52
knm$edge.par[[43]][1,1,1] <- 53
knm$edge.par[[43]][2,2,1] <- 53
knm$edge.par[[44]][1,1,1] <- 54
knm$edge.par[[44]][2,2,1] <- 54
knm$edge.par[[45]][1,1,1] <- 55
knm$edge.par[[45]][2,2,1] <- 55

#set.seed(6)
knm$par <- runif(55,-1.5,1.1)
knm$par # "true" theta
out.pot <- make.pots(parms = knm$par,  crf = knm,  rescaleQ = T, replaceQ = T)

# So now sample from the model as if we obtained an experimental sample:
num.samps <- 1000
#set.seed(1)
samps <- sample.exact(knm, num.samps)
mrf.sample.plot(samps)


psl <- make.crf(adj, n.states)
psl <- make.features(psl)
psl <- make.par(psl, 55)
psl$n.edges

psl$node.par[1,1,] <- 1
psl$node.par[2,1,] <- 2
psl$node.par[3,1,] <- 3
psl$node.par[4,1,] <- 4
psl$node.par[5,1,] <- 5
psl$node.par[6,1,] <- 6
psl$node.par[7,1,] <- 7
psl$node.par[8,1,] <- 8
psl$node.par[9,1,] <- 9
psl$node.par[10,1,] <- 10
psl$edge.par[[1]][1,1,1] <- 11
psl$edge.par[[1]][2,2,1] <- 11
psl$edge.par[[2]][1,1,1] <- 12
psl$edge.par[[2]][2,2,1] <- 12
psl$edge.par[[3]][1,1,1] <- 13
psl$edge.par[[3]][2,2,1] <- 13
psl$edge.par[[4]][1,1,1] <- 14
psl$edge.par[[4]][2,2,1] <- 14
psl$edge.par[[5]][1,1,1] <- 15
psl$edge.par[[5]][2,2,1] <- 15
psl$edge.par[[6]][1,1,1] <- 16
psl$edge.par[[6]][2,2,1] <- 16
psl$edge.par[[7]][1,1,1] <- 17
psl$edge.par[[7]][2,2,1] <- 17
psl$edge.par[[8]][1,1,1] <- 18
psl$edge.par[[8]][2,2,1] <- 18
psl$edge.par[[9]][1,1,1] <- 19
psl$edge.par[[9]][2,2,1] <- 19
psl$edge.par[[10]][1,1,1] <- 20
psl$edge.par[[10]][2,2,1] <- 20

psl$edge.par[[11]][1,1,1] <- 21
psl$edge.par[[11]][2,2,1] <- 21
psl$edge.par[[12]][1,1,1] <- 22
psl$edge.par[[12]][2,2,1] <- 22
psl$edge.par[[13]][1,1,1] <- 23
psl$edge.par[[13]][2,2,1] <- 23
psl$edge.par[[14]][1,1,1] <- 24
psl$edge.par[[14]][2,2,1] <- 24
psl$edge.par[[15]][1,1,1] <- 25
psl$edge.par[[15]][2,2,1] <- 25
psl$edge.par[[16]][1,1,1] <- 26
psl$edge.par[[16]][2,2,1] <- 26
psl$edge.par[[17]][1,1,1] <- 27
psl$edge.par[[17]][2,2,1] <- 27
psl$edge.par[[18]][1,1,1] <- 28
psl$edge.par[[18]][2,2,1] <- 28
psl$edge.par[[19]][1,1,1] <- 29
psl$edge.par[[19]][2,2,1] <- 29
psl$edge.par[[20]][1,1,1] <- 30
psl$edge.par[[20]][2,2,1] <- 30

psl$edge.par[[21]][1,1,1] <- 31
psl$edge.par[[21]][2,2,1] <- 31
psl$edge.par[[22]][1,1,1] <- 32
psl$edge.par[[22]][2,2,1] <- 32
psl$edge.par[[23]][1,1,1] <- 33
psl$edge.par[[23]][2,2,1] <- 33
psl$edge.par[[24]][1,1,1] <- 34
psl$edge.par[[24]][2,2,1] <- 34
psl$edge.par[[25]][1,1,1] <- 35
psl$edge.par[[25]][2,2,1] <- 35
psl$edge.par[[26]][1,1,1] <- 36
psl$edge.par[[26]][2,2,1] <- 36
psl$edge.par[[27]][1,1,1] <- 37
psl$edge.par[[27]][2,2,1] <- 37
psl$edge.par[[28]][1,1,1] <- 38
psl$edge.par[[28]][2,2,1] <- 38
psl$edge.par[[29]][1,1,1] <- 39
psl$edge.par[[29]][2,2,1] <- 39
psl$edge.par[[30]][1,1,1] <- 40
psl$edge.par[[30]][2,2,1] <- 40

psl$edge.par[[31]][1,1,1] <- 41
psl$edge.par[[31]][2,2,1] <- 41
psl$edge.par[[32]][1,1,1] <- 42
psl$edge.par[[32]][2,2,1] <- 42
psl$edge.par[[33]][1,1,1] <- 43
psl$edge.par[[33]][2,2,1] <- 43
psl$edge.par[[34]][1,1,1] <- 44
psl$edge.par[[34]][2,2,1] <- 44
psl$edge.par[[35]][1,1,1] <- 45
psl$edge.par[[35]][2,2,1] <- 45
psl$edge.par[[36]][1,1,1] <- 46
psl$edge.par[[36]][2,2,1] <- 46
psl$edge.par[[37]][1,1,1] <- 47
psl$edge.par[[37]][2,2,1] <- 47
psl$edge.par[[38]][1,1,1] <- 48
psl$edge.par[[38]][2,2,1] <- 48
psl$edge.par[[39]][1,1,1] <- 49
psl$edge.par[[39]][2,2,1] <- 49
psl$edge.par[[40]][1,1,1] <- 50
psl$edge.par[[40]][2,2,1] <- 50

psl$edge.par[[41]][1,1,1] <- 51
psl$edge.par[[41]][2,2,1] <- 51
psl$edge.par[[42]][1,1,1] <- 52
psl$edge.par[[42]][2,2,1] <- 52
psl$edge.par[[43]][1,1,1] <- 53
psl$edge.par[[43]][2,2,1] <- 53
psl$edge.par[[44]][1,1,1] <- 54
psl$edge.par[[44]][2,2,1] <- 54
psl$edge.par[[45]][1,1,1] <- 55
psl$edge.par[[45]][2,2,1] <- 55

Delta.alpha.info <- delta.alpha(crf = psl, samples = samps, printQ = F)
Delta.alpha <- Delta.alpha.info$Delta.alpha
Delta.alpha

y <-c(samps[,1], samps[,2], samps[,3], samps[,4], samps[,5], samps[,6], samps[,7], samps[,8], samps[,9], samps[,10])
y
y[which(y==2)] <- 0
y
length(y)
dim(Delta.alpha)

Ma <- glm(y ~ Delta.alpha[,1] +
              Delta.alpha[,2] +
              Delta.alpha[,3] +
              Delta.alpha[,4] +
              Delta.alpha[,5] +
              Delta.alpha[,6] +
              Delta.alpha[,7] +
              Delta.alpha[,8] +
              Delta.alpha[,9] +
              Delta.alpha[,10] +
            Delta.alpha[,11] +
            Delta.alpha[,12] +
            Delta.alpha[,13] +
            Delta.alpha[,14] +
            Delta.alpha[,15] +
            Delta.alpha[,16] +
            Delta.alpha[,17] +
            Delta.alpha[,18] +
            Delta.alpha[,19] +
            Delta.alpha[,20] +
            Delta.alpha[,21] +
            Delta.alpha[,22] +
            Delta.alpha[,23] +
            Delta.alpha[,24] +
            Delta.alpha[,25] +
            Delta.alpha[,26] +
            Delta.alpha[,27] +
            Delta.alpha[,28] +
            Delta.alpha[,29] +
            Delta.alpha[,30] +
            Delta.alpha[,31] +
            Delta.alpha[,32] +
            Delta.alpha[,33] +
            Delta.alpha[,34] +
            Delta.alpha[,35] +
            Delta.alpha[,36] +
            Delta.alpha[,37] +
            Delta.alpha[,38] +
            Delta.alpha[,39] +
            Delta.alpha[,40] +
            Delta.alpha[,41] +
            Delta.alpha[,42] +
            Delta.alpha[,43] +
            Delta.alpha[,44] +
            Delta.alpha[,45] +
            Delta.alpha[,46] +
            Delta.alpha[,47] +
            Delta.alpha[,48] +
            Delta.alpha[,49] +
            Delta.alpha[,50] +
            Delta.alpha[,51] +
            Delta.alpha[,52] +
            Delta.alpha[,53] +
            Delta.alpha[,54] +
            Delta.alpha[,55] - 1,
          family=binomial(link="logit"))
summary(Ma)

cbind(coef(Ma), knm$par)

psl$par <- as.numeric(coef(Ma))

out.pot2 <- make.pots(parms = psl$par,  crf = psl,  rescaleQ = T, replaceQ = T)
psl$node.pot
psl$edge.pot

# Node and edge beliefs:
psl.bel <- infer.exact(psl)
knm.bel <- infer.exact(knm)

cbind(psl.bel$node.bel[,1], knm.bel$node.bel[,1])
cbind(psl.bel$node.bel[,2], knm.bel$node.bel[,2])

psl.bel$edge.bel[[1]]
knm.bel$edge.bel[[1]]
psl.bel$edge.bel[[2]]
knm.bel$edge.bel[[2]]
psl.bel$edge.bel[[3]]
knm.bel$edge.bel[[3]]
psl.bel$edge.bel[[4]]
knm.bel$edge.bel[[4]]
psl.bel$edge.bel[[5]]
knm.bel$edge.bel[[5]]
psl.bel$edge.bel[[6]]
knm.bel$edge.bel[[6]]
psl.bel$edge.bel[[7]]
knm.bel$edge.bel[[7]]
psl.bel$edge.bel[[8]]
knm.bel$edge.bel[[8]]
psl.bel$edge.bel[[9]]
knm.bel$edge.bel[[9]]
psl.bel$edge.bel[[10]]
knm.bel$edge.bel[[10]]

psl.bel$edge.bel[[21]]
knm.bel$edge.bel[[21]]
psl.bel$edge.bel[[22]]
knm.bel$edge.bel[[22]]
psl.bel$edge.bel[[23]]
knm.bel$edge.bel[[23]]
psl.bel$edge.bel[[24]]
knm.bel$edge.bel[[24]]
psl.bel$edge.bel[[25]]
knm.bel$edge.bel[[25]]
psl.bel$edge.bel[[26]]
knm.bel$edge.bel[[26]]
psl.bel$edge.bel[[27]]
knm.bel$edge.bel[[27]]
psl.bel$edge.bel[[28]]
knm.bel$edge.bel[[28]]
psl.bel$edge.bel[[29]]
knm.bel$edge.bel[[29]]
psl.bel$edge.bel[[30]]
knm.bel$edge.bel[[30]]

psl.bel$edge.bel[[41]]
knm.bel$edge.bel[[41]]
psl.bel$edge.bel[[42]]
knm.bel$edge.bel[[42]]
psl.bel$edge.bel[[43]]
knm.bel$edge.bel[[43]]
psl.bel$edge.bel[[44]]
knm.bel$edge.bel[[44]]
psl.bel$edge.bel[[45]]
knm.bel$edge.bel[[45]]

# Configuration probabilities:
options("max.print" = 20000)

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

colnames(joint.dist.info)[c(9,10,8:1)]


psl.cp <- round(joint.dist.info[,11]*100, 5)
knm.cp <- round(joint.dist.info.knm[,11]*100, 5)
cbind(joint.dist.info[,c(9,10,8:1)], psl.cp, knm.cp)

sum(psl.cp)
sum(knm.cp)

plot(psl.cp, typ="h")
plot(knm.cp, typ="h")

