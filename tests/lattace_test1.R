jmat <- t(array(c("A","B","C","D","E","FF","G","H","II","J","K","L","M","N","O","P") ,c(4,4)))
jmat

eq <- NULL
# Horizontal links
for(i in 1:nrow(jmat)){
  for(j in 1:(ncol(jmat)-1)){
    #print(i)
    htrm <- paste0(jmat[i,j], ":", jmat[i,j+1])
    print(htrm)
    eq <- c(eq, htrm)
  }
}

# Vertical links
for(j in 1:ncol(jmat)){
  for(i in 1:(nrow(jmat)-1)) {
    vtrm <- paste0(jmat[i,j], ":", jmat[i+1,j])
    print(vtrm)
    eq <- c(eq, vtrm)
  }
}

#Cross links
for(i in 1:(nrow(jmat)-1)) {
  for(j in 1:ncol(jmat)){
    if((j != 1)&(j != ncol(jmat))){
      #print(paste(jmat[i,j], "is NOT a corner"))
      #print(paste("   ", jmat[i+1,j-1], "is a neighbor"))
      #print(paste("   ", jmat[i+1,j],   "is a neighbor"))
      #print(paste("   ", jmat[i+1,j+1], "is a neighbor"))

      cterm1 <- paste0(jmat[i,j], ":", jmat[i+1,j-1])
      cterm2 <- paste0(jmat[i,j], ":", jmat[i+1,j])
      cterm3 <- paste0(jmat[i,j], ":", jmat[i+1,j+1])
      print(cterm1)
      print(cterm2)
      print(cterm3)
      eq <- c(eq, cterm1, cterm2, cterm3)
    } else {
      #print(paste(jmat[i,j], "is a corner"))
      if(j == 1) {
        #print("    Start corner")
        #print(paste("   ", jmat[i+1,j],   "is a neighbor"))
        #print(paste("   ", jmat[i+1,j+1], "is a neighbor"))

        cterm2 <- paste0(jmat[i,j], ":", jmat[i+1,j])
        cterm3 <- paste0(jmat[i,j], ":", jmat[i+1,j+1])
        print(cterm2)
        print(cterm3)
        eq <- c(eq, cterm2, cterm3)
      } else {
        #print("    Stop corner")
        #print(paste("   ", jmat[i+1,j-1], "is a neighbor"))
        #print(paste("   ", jmat[i+1,j],   "is a neighbor"))

        cterm1 <- paste0(jmat[i,j], ":", jmat[i+1,j-1])
        cterm2 <- paste0(jmat[i,j], ":", jmat[i+1,j])
        print(cterm1)
        print(cterm2)
        eq <- c(eq, cterm1, cterm2)
      }
    }
  }
  print("---------")
}

eq

#unique(eq)

grphf <- paste0("~", eq[1], sep="")
for(i in 2:length(eq)){
  grphf <- paste0(grphf, "+", eq[i], sep="")
}
grphf <- as.formula(grphf)


adj        <- ug(grphf, result="matrix") # adjacency (connection) matrix
node.names <- colnames(adj)
# Check the graph:
gp <- ug(grphf, result = "graph")
plot(gp)
dev.off()
