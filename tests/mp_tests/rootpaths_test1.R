library(CRFutil)

# Model:
grphf <- ~1:3 + 2:3 + 3:4 + 3:5 + 5:6 + 6:7
gamt <-  rbind(
  c(1,3),
  c(2,3),
  c(3,4),
  c(3,5),
  c(5,6),
  c(6,7))

adj <- edges2adj(gamt, plotQ = T)

knm <- make.empty.field(
  graph.eq             = NULL,
  adj.mat              = adj,
  parameterization.typ = "ising2",
  node.par             = NULL,
  edge.par             = NULL,
  plotQ                = T)
dump.crf(knm)
knm$adj.nodes[[5]] # whle not bottom??
knm$adj.nodes[[3]] # Call explore at each level and return ne nodes. If node is a leaf, stop
knm$adj.nodes[[6]]

#Call explore
#if leaf node stop
# else loop over neighbor list


ptht <- sp.between(ug(grphf), "5", "4")
ptht[[1]]$path_detail

# Leaf nodes:
which(sapply(1:knm$n.nodes,function(xx){length(knm$adj.nodes[[xx]])}) == 1)


# Find all paths from root to leaf nodes
ugp <- ug(grphf)
leaf.nodes <- leaves(ugp)
rt <- 4 # root node
#sp.between(ugp, as.character(rt), as.character(leaf.nodes[1]))[[1]]$path_detail
lapply(1:length(leaf.nodes), function(xx){sp.between(ugp, as.character(rt), leaf.nodes[xx])[[1]]$path_detail})
# If leaf node is root, it will be its own path. drop it somehow

sp.between(ugp, as.character(rt), as.character(leaf.nodes[1]))[[1]]

#library(RBGL)
get.path(ugp, 4, 7)
get.root.paths(grphf,root.node = 5)
