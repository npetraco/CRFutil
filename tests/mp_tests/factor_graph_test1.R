library(CRFutil)

# Model:
grphf <- ~1:3 + 2:3 + 3:4 + 3:5 + 5:6 + 6:7

ugm <- ug(grphf)
plot(ugm)

# Root paths to leaves:
get.root.paths(grphf,root.node = 6)

# Build pair-wise factor graph from UG model
adj <- as(ugm, "matrix")
adj.uptri <- upper.tri(adj) * adj
edg.mat <- which(adj.uptri == 1, arr.ind = T)
edg.mat <- cbind(
  as.numeric(rownames(adj)[edg.mat[,1]]),
  as.numeric(colnames(adj)[edg.mat[,2]])
)
edg.mat <- t(sapply(1:nrow(edg.mat), function(xx){sort(edg.mat[xx,])}))
edg.mat

#edg.mat2 <- array(0,c(2*nrow(edg.mat), 2))
edg.mat2 <- NULL
# edge factors
for(i in 1:nrow(edg.mat)) {
  factor.nme <- paste0("f",edg.mat[i,1],edg.mat[i,2])
  edg.mat2 <- rbind(
    edg.mat2,
    c(edg.mat[i,1], factor.nme),
    c(factor.nme, edg.mat[i,2])
  )
}
edg.mat2

# tack on node factors
edg.mat2 <- rbind(
  cbind(
    paste0("f",sort(ugm@nodes)),
    sort(ugm@nodes)
  ),
  edg.mat2
)
edg.mat2


library(igraph)
g <- graph_from_data_frame(data.frame(edg.mat2), directed = FALSE)
#plot(g)
nde.nms <- V(g)$name
V(g)$type <- sapply(1:length(nde.nms), function(xx){length(strsplit(nde.nms[xx], split = "f")[[1]])})

# Plot factor graph:
cols <- c("steelblue", "red")
shps <- c("circle", "square")
plot(g,
     vertex.color = cols[as.numeric(V(g)$type)],
     vertex.shape = shps[as.numeric(V(g)$type)]
)

g2 <- igraph.to.graphNEL(g)
plot(g2)

rpths <- get.root.paths(g2,root.node = "5")
sched.mat <- array(NA,c(length(rpths), max(sapply(1:length(rpths), function(xx){length(rpths[[xx]])}))))
for(i in 1:length(rpths)) {
  rev.pth <- rev(rpths[[i]])
  for(j in 1:length(rpths[[i]])) {
    #print(rpths[[i]][j])
    sched.mat[i,j] <- rev.pth[j]
  }
}
sched.mat
t(sched.mat)
