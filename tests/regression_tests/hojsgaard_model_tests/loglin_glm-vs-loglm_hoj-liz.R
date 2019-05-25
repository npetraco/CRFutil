library(gRbase)
library(Rgraphviz)
library(MASS)
library(gRim)

# Observed configurations: X
data("lizardRAW")
data("lizardAGG")
X  <- lizardRAW   # Observed states
Xf <- lizardAGG   # Freq table of observed states

# Raw case list to Aggregated Frequency Table
as.data.frame(ftable(lizardRAW))
lizardAGG

## Raw case-list to contingency table (Observed states)
xtabs(~., data=lizardRAW)


# contingency table
lizard <- xtabs(~., data=lizardRAW)
dim(lizard)
Xc <- lizard     # Fold into a contingency table

# in gRim:
grphf <- ~species:diam + species:height
plot(ug(grphf))

ll3  <- dmod(grphf, data = X)
class(ll3)
ll3

# using loglm
ll2 <- loglm(Freq ~ species:diam + species:height, data = Xf)
summary(ll2)
coef(ll2)

# Model matrix ??????
Xf
Xf[,1]
dia.dc <- Xf[,1]
contrasts(dia.dc)
contrasts(dia.dc) <- contr.sum(2)
contrasts(dia.dc)

Xf[,2]
hgt.dc <- Xf[,2]
contrasts(hgt.dc)
contrasts(hgt.dc) <- contr.sum(2)
contrasts(hgt.dc)

Xf[,3]
spc.dc <- Xf[,3]
contrasts(spc.dc)
contrasts(spc.dc) <- contr.sum(2)
contrasts(spc.dc)

Xm <- model.matrix(~spc.dc + dia.dc + hgt.dc + spc.dc:dia.dc + spc.dc:hgt.dc,
             contrasts = list(spc.dc = "contr.sum",
                              dia.dc = "contr.sum",
                              hgt.dc = "contr.sum"))

loglm( ~ species:diam + species:height, data = Xc)
loglm( ~ species + diam + height + species:diam + species:height, data = Xc)
loglm(Freq~spc.dc + dia.dc + hgt.dc + spc.dc:dia.dc + spc.dc:hgt.dc, data = Xf)
#loglm( ~ species:diam + species:height, data = X)
#loglm(Xf[,4] ~ Xm[,2] + Xm[,3] + Xm[,4] + Xm[,5] + Xm[,6])
#loglm(Freq~Xm[,2:6], data = Xf)

c1 <- coef(loglm( ~ species:diam + species:height, data = Xc))
c2 <- coef(loglm( ~ species + diam + height + species:diam + species:height, data = Xc))
c3 <- coef(loglm(Freq~spc.dc + dia.dc + hgt.dc + spc.dc:dia.dc + spc.dc:hgt.dc, data = Xf))

c1[[1]]
c2[[1]]
c3[[1]]

c1[[2]]
c2[[2]]
c3[[3]] # ??

c1[[3]]
c2[[3]]
c3[[4]] # ??

c1[[4]]
c2[[4]]
c3[[2]] # ??

c1[[5]]
c2[[5]]
t(c3[[5]])

c1[[6]]
c2[[6]]
t(c3[[6]])

# glm next
glm.f <- glm(Freq~spc.dc + dia.dc + hgt.dc + spc.dc:dia.dc + spc.dc:hgt.dc, data = Xf, family = poisson(link="log"))
names(c1)
lmcoefs <- c(c1[[1]][1], c1[[2]][1], c1[[3]][1], c1[[4]][1], c1[[5]][1,1], c1[[6]][1,1]) # loglm coefs
names(lmcoefs) <- names(c1)

lmcoefs     # compare
coef(glm.f)

# Default contrast coding: treatment coding
glm.reg <- glm(Freq ~ species + diam + height + species:diam + species:height, data = Xf, family = poisson(link="log"))
summary(glm.reg)
coef(glm.reg) - coef(glm.reg)[1] # Nope...

# Model matrix for treatment coding:
Xmdf <- model.matrix(~species + diam + height + species:diam + species:height, data = Xf)
Xmdf

?contrasts
#Xf[,3]
#spc.dcf <- Xf[,3]
colnames(Xmdf)
glm(Xf[,4] ~ Xmdf[,2] + Xmdf[,3] + Xmdf[,4] + Xmdf[,5] + Xmdf[,6], family = poisson(link="log"))
coef(glm.reg)

Xmdf
Xm

# This works. craft model around X also????
# Would it be equivalent to a logistic for the whole state?
model.matrix(~species + diam + height + species:diam + species:height, data = X)
