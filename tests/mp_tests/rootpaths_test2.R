library(CRFutil)

# Model:
grphf <- ~1:3 + 2:3 + 3:4 + 3:5 + 5:6 + 6:7

# Find all paths from root to leaves:
get.root.paths(grphf,root.node = 5)


# Model with a loop:
grphf <- ~1:2 + 2:3 + 3:4 + 4:1 + 1:5
plot(ug(grphf))
get.root.paths(grphf,root.node = 2)
