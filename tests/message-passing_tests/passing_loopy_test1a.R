library(CRFutil)
library(Rgraphviz)

# Loopy Message pass test

# Model:
grf <- ~A:B
adj <- ug(grf, result="matrix")
f0  <- function(y){ as.numeric(c((y==1),(y==2)))}

n.states <- 2
km       <- make.crf(adj, n.states)

# Node pots:
PsiA <- c(3,1)
PsiB <- c(3,1)

# Edge pots:
PsiAB <-
  rbind(
    c(3, 6.1),
    c(6.1, 3)
  )

km$node.pot[1,]  <- PsiA
km$node.pot[2,]  <- PsiB
km$edge.pot[[1]] <- PsiAB

# Dress potentials for message passing:
gr.pots  <- make.gRbase.potentials(crf = km, node.names = rownames(adj), state.nmes = c("up", "dn"))
ep       <- gr.pots$edge.potentials # Edge potentials in gRbase format
np       <- gr.pots$node.potentials # Node potentials in gRbase format
all.pots <- union(np, ep)           # Put all node and edge potentials into one list so we don't have to search though two lists for them

nd.nms          <- nodes(ug(grf))                                               # node names
edg.nms         <- cbind(nd.nms[km$edges[,1]], nd.nms[km$edges[,2]])            # edge names
nd.pot.nms      <- paste0("f", nd.nms)                                          # node potential names
edg.pot.nms     <- paste0("f", nd.nms[km$edges[,1]], "-", nd.nms[km$edges[,2]]) # edge potential names
names(all.pots) <- c(nd.pot.nms, edg.pot.nms) # Put the same names on the potentials as will appear in the schedule
all.pots        # Potentials dressed up in gRbase format and ready for message passing


# Convert MRF to a pair-wise factor graph
pwfg <- mrf2pwfg(grf, plotQ=T)
dev.off()
plot(pwfg, nodeAttrs=makeNodeAttrs(pwfg, fontsize=30))

# Try a loopy passes??
msg.cont.info <- init.loopy.message.storage(pwfg)
msg.bxs       <- msg.cont.info$message.container
msg.sch       <- msg.cont.info$message.schedule.mat

msg.bxs # message storage list (mailboxes)
msg.sch # schedule


# Loopy Pass messages according to schedule:
max.iter <- 2
for(iter in 1:max.iter){

  for(i in 1:nrow(msg.sch)){
    msg.num <- i
    print(msg.sch[msg.num,]) # message info
    msg.nd.nmes <- as.character(msg.sch[msg.num,c(2,3)]) # start-end nodes of message

    msg.info <- pass.message(start.node = msg.nd.nmes[1],
                             end.node   = msg.nd.nmes[2],
                             factor.graph = pwfg,
                             pots.list  = all.pots,
                             mailboxes.list = msg.bxs,
                             printQ = F)

    #print(msg.info)

    msg.bxs[[msg.info$mailbox.idx]] <- msg.info$message

  }

  print(paste("===== Loopy iteration:", iter, "done ====="))

}
msg.bxs


# Get exact marginals to check against:
all.marginals <- infer.exact(km)
all.marginals

# Node marginals
nde.nms <- nodes(ug(grf))
for(i in 1:length(nde.nms)){
  ndm <- node.marginal(nde.nms[i], pwfg, msg.bxs)
  print(paste0("bel(",nde.nms[i],"):"))
  print(ndm)
  print(round(ndm - all.marginals$node.bel[i,], 3)) # Diffs = about 0??
  print("==================================")
}


# Edge marginals
for(i in 1:nrow(km$edges)){
  lnd <- nd.nms[km$edges[i,1]]
  rnd <- nd.nms[km$edges[i,2]]
  print(paste0("bel(",lnd,",",rnd,"):"))
  edgm <- edge.marginal(v.start.node = lnd, v.end.node = rnd,
                        factor.graph = pwfg, pots.list = all.pots, mailbox.list = msg.bxs)
  print(edgm)
  print(round(edgm - all.marginals$edge.bel[[i]], 3)) # Diffs = about 0??
  print("==================================")
}
