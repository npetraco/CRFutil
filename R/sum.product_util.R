#' Get a (shortest) path between two nodes
#'
#' Get a (shortest) path between two nodes
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
get.path <- function(a.graphnel.obj, start.node, end.node){ # **** BROKEN?? SOME SEGMENTS MISSED IN SOME POLYTREES ******

  path.info <- sp.between(a.graphnel.obj, as.character(start.node), as.character(end.node))[[1]]
  if(is.na(path.info[1])) {
    #print("No path!")
    sepath <- NULL
  } else {
    sepath <- path.info$path_detail
  }

  return(sepath)
}


#' Get all paths from chosen root node to all leaf nodes
#'
#' Get all paths from chosen root node to all leaf nodes
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
get.root.paths <- function(graph.obj, root.node, serial.schedsQ=FALSE) {

  if(class(graph.obj) == "formula") {
    und.gph    <- ug(graph.obj)
  } else if(class(graph.obj) == "graphNEL") {
    und.gph    <- graph.obj
  } else {
    stop("Must enter a formula or graphNEL object for the graph.obj arguement.")
  }
  leaf.nodes <- leaves(und.gph)

  all.root.pths <- lapply(1:length(leaf.nodes), function(xx){get.path(und.gph, root.node, leaf.nodes[xx])})

  # Check and see if any of the paths are NULL. Probably chose a leaf node as root in that case:
  null.pthsQ <- sapply(1:length(all.root.pths), function(xx){is.null(all.root.pths[[xx]])})
  drp.pths   <- which(null.pthsQ == TRUE)
  if(length(drp.pths) !=0) {
    all.root.pths <- all.root.pths[-drp.pths]
  }

  if(serial.schedsQ == TRUE) {
    sched.mat.forward  <- array(NA,c(length(all.root.pths), max(sapply(1:length(all.root.pths), function(xx){length(all.root.pths[[xx]])}))))
    sched.mat.backward <- array(NA,c(length(all.root.pths), max(sapply(1:length(all.root.pths), function(xx){length(all.root.pths[[xx]])}))))
    for(i in 1:length(all.root.pths)) {
      rev.pth <- rev(all.root.pths[[i]])
      pth     <- all.root.pths[[i]]
      for(j in 1:length(all.root.pths[[i]])) {
        sched.mat.forward[i,j]  <- rev.pth[j]
        sched.mat.backward[i,j] <- pth[j]
      }
    }
    all.root.pths <- list(sched.mat.forward, sched.mat.backward)
    names(all.root.pths) <- c("forward", "backward")
  }

  return(all.root.pths)

}

#' Convert Markov Random Field to pairwise (two-body) factor graph:
#'
#' Convert Markov Random Field to pairwise (two-body) factor graph:
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
mrf2pwfg <- function(graph.obj, plotQ=FALSE) {  # BROKEN!!!!!! LARGE POLYTREES WITH MANY NODES MISCONNECTING NODES

  if(class(graph.obj) == "formula") {
    und.gph    <- ug(graph.obj)
  } else if(class(graph.obj) == "graphNEL") {
    und.gph    <- graph.obj
  } else {
    stop("Must enter a formula or graphNEL object for the graph.obj arguement.")
  }

  # Use igraph functionality to build pair-wise factor graph from UG model
  adjmat <- as(und.gph, "matrix")
  adjmat.uptri <- upper.tri(adjmat) * adjmat

  # edges:
  edgmat <- which(adjmat.uptri == 1, arr.ind = T)
  edgmat <- cbind(
    as.numeric(rownames(adjmat)[edgmat[,1]]),  # This means that input graph nodes must be numbers....
    as.numeric(colnames(adjmat)[edgmat[,2]])   # Can we remove this restriction or try to get around it?
  )
  edgmat <- t(sapply(1:nrow(edgmat), function(xx){sort(edgmat[xx,])}))

  edgmat2 <- NULL
  # add edge factors between connected nodes:
  for(i in 1:nrow(edgmat)) {
    factor.nme <- paste0("f",edgmat[i,1],edgmat[i,2])  # PROBLEM HERE? EG f12 edge name gets confused with f12 node name
    edgmat2 <- rbind(
      edgmat2,
      c(edgmat[i,1], factor.nme),
      c(factor.nme, edgmat[i,2])
    )
  }

  # tack on node factors:
  edgmat2 <- rbind(
    cbind(
      paste0("f",sort(und.gph@nodes)),
      sort(und.gph@nodes)
    ),
    edgmat2
  )

  fg <- graph_from_data_frame(data.frame(edgmat2), directed = FALSE)
  nde.nms <- V(fg)$name
  V(fg)$type <- sapply(1:length(nde.nms), function(xx){length(strsplit(nde.nms[xx], split = "f")[[1]])})

  if(plotQ==TRUE){
    # Plot factor graph:
    cols <- c("steelblue", "red")
    shps <- c("circle", "square")
    plot(fg,
         vertex.color = cols[as.numeric(V(fg)$type)],
         vertex.shape = shps[as.numeric(V(fg)$type)]
    )
  }

  # Convert back to graphNEL format:
  fg <- igraph.to.graphNEL(fg)
  return(fg)
}


#' Convert Markov Random Field to pairwise (two-body) factor graph:
#'
#' Convert Markov Random Field to pairwise (two-body) factor graph:
#'
#' This version allows for nodes names to be characters. It just spits out a warning if they are.
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
mrf2pwfg2 <- function(graph.obj, plotQ=FALSE) {  # mrf2pwfg BROKEN!!!!!! Re-write try


  if(class(graph.obj) == "formula") {
    und.gph    <- ug(graph.obj)
  } else if(class(graph.obj) == "graphNEL") {
    und.gph    <- graph.obj
  } else {
    stop("Must enter a formula or graphNEL object for the graph.obj arguement.")
  }

  # Node names may not be numbers. For this function, we really need the order-ability of numbers.
  # So if node names are not numbers, assign each name a number
  node.names          <- und.gph@nodes
  node.names.numericQ <- !(NA %in% as.numeric(node.names))
  orig.node.names     <- node.names
  if(node.names.numericQ == FALSE){

    warning("Node names are not numeric. Assigning nodes a cononical order.")

    und.gph <- graph_from_graphnel(und.gph)
    V(und.gph)$name <- 1:length(node.names)
    und.gph <- as_graphnel(und.gph)
  }

  # Use igraph functionality to build pair-wise factor graph from UG model
  adjmat <- as(und.gph, "matrix")
  adjmat.uptri <- upper.tri(adjmat) * adjmat

  # edges:
  edgmat <- which(adjmat.uptri == 1, arr.ind = T)
  edgmat <- cbind(
    as.numeric(rownames(adjmat)[edgmat[,1]]),  # This means that input graph nodes must be numbers....
    as.numeric(colnames(adjmat)[edgmat[,2]])   # Can we remove this restriction or try to get around it?
  )
  edgmat <- t(sapply(1:nrow(edgmat), function(xx){sort(edgmat[xx,])}))

  # In case orig node names were characters
  edgmatC <- cbind(orig.node.names[edgmat[,1]], orig.node.names[edgmat[,2]])

  edgmat2  <- NULL
  edgmat2C <- NULL # In case orig node names were characters
  # add edge factors between connected nodes:
  for(i in 1:nrow(edgmat)) {
    factor.nme  <- paste0("f",edgmat[i,1],"-",edgmat[i,2])   # PROBLEM WAS HERE. EG f12 edge name gets confused with f12 node name. Added in the - to the name
    factor.nmeC <- paste0("f",edgmatC[i,1],"-",edgmatC[i,2]) # In case orig node names were characters

    edgmat2 <- rbind(
      edgmat2,
      c(edgmat[i,1], factor.nme),
      c(factor.nme, edgmat[i,2])
    )

    # In case orig node names were characters
    edgmat2C <- rbind(
      edgmat2C,
      c(edgmatC[i,1], factor.nmeC),
      c(factor.nmeC, edgmatC[i,2])
    )

  }

  # tack on node factors:
  edgmat2 <- rbind(
    cbind(
      paste0("f",sort(und.gph@nodes)),
      sort(und.gph@nodes)
    ),
    edgmat2
  )
  # In case orig node names were characters
  edgmat2C <- rbind(
    cbind(
      paste0("f",orig.node.names[order(und.gph@nodes)]),
      orig.node.names[order(und.gph@nodes)]
    ),
    edgmat2C
  )


  fg <- graph_from_data_frame(data.frame(edgmat2), directed = FALSE)
  # In case orig node names were characters
  fgC <- graph_from_data_frame(data.frame(edgmat2C), directed = FALSE)
  print(cbind(edgmat2, edgmat2C))

  nde.nms     <- V(fg)$name
  V(fg)$type  <- sapply(1:length(nde.nms), function(xx){length(strsplit(nde.nms[xx], split = "f")[[1]])})
  node.types  <- V(fg)$type # For plotting in case node names are converted to characters
  V(fgC)$type <- node.types # In case orig node names were characters

  if(plotQ==TRUE){
    # Plot factor graph:
    cols <- c("steelblue", "red")
    shps <- c("circle", "rectangle")
    plot(fgC,
         vertex.color = cols[as.numeric(node.types)],
         vertex.shape = shps[as.numeric(node.types)]
    )
  }

  # Convert back to graphNEL format:
  fg  <- igraph.to.graphNEL(fg)  # Should we return the ordinal node name version too?? We had been doing that....
  fgC <- igraph.to.graphNEL(fgC) # In case orig node names were characters

  return(fgC)

}

#' Make factor to variable (node) message
#'
#' Make factor to variable (node) message
#' \mu_{f\rightarrow X} = \sum_{\backslash X} \Big( f(X) \prod_{Y \in \text{ne}(f)\backslash X} \mu_{Y\rightarrow f} \Big)
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
make.f2v.msg <- function(in.v.msgs.list = NULL, f.msg, out.v.nme){

  if(is.null(in.v.msgs.list)) {
    # v.nme check for initalization??
    msg.out <- f.msg # For initialization case
  } else {
    msg.prod <- f.msg %a*% ar_prod_list(in.v.msgs.list)
    msg.out <- ar_marg(msg.prod, out.v.nme)
  }

  return(msg.out) # Guarantee sending this out as a list?????

}


#' Make factor to variable (node) message
#'
#' Make variable (node) to factormessage
#' XXXXXX
#'
#' The function will XXXX
#'
#' @param XX The XX
#' @return The function will XX
#'
#'
#' @export
make.v2f.msg <- function(in.f.msgs.list = NULL){

  if(is.null(in.f.msgs.list)) {
    msg.out <- NULL # For initialization case. Was "id".
  } else {
    msg.out <- ar_prod_list(in.f.msgs.list)
  }

  return(msg.out)  # Guarantee sending this out as a list?????

}


#' Initialize an empty set of list containers to hold messages
#'
#' Initialize an empty set of list containers to hold messages
#' A NULL message can be interpreted as an id message if necessary.
#'
#' The function will XXXX
#'
#' @param a.schedule A message passing schedule to a root node.
#' @return The function will XX
#'
#'
#' @export
init.message.storage <- function(a.schedule){

  a.mailroom        <- NULL
  message.names     <- NULL
  message.names.mat <- NULL

  # Determine message passing schedule first:

  # Initialize forward messages:
  for(i in 1:nrow(a.schedule$forward)){
    for(j in 1:(ncol(a.schedule$forward) - 1)){

      if(is.na(a.schedule$forward[i, j+1])) {
        break()
      } else {
        msg.symb      <- paste0(a.schedule$forward[i, j], ".", a.schedule$forward[i, j+1])
        message.names <- c(message.names, msg.symb)

        # Another format for the messages passed. Includes the schedule index j
        message.names.mat <- rbind(message.names.mat, c(j, a.schedule$forward[i, j], a.schedule$forward[i, j+1], msg.symb))
      }

    }
  }

  # Initialize backward messages:
  back.pass.off <- ncol(a.schedule$forward) # Offset index for backwards pass schedule
  for(i in 1:nrow(a.schedule$backward)){
    for(j in 1:(ncol(a.schedule$backward) - 1)){

      if(is.na(a.schedule$backward[i, j+1])) {
        break()
      } else {

        msg.symb      <- paste0(a.schedule$backward[i, j], ".", a.schedule$backward[i, j+1])
        message.names <- c(message.names, msg.symb)

        # Another format for the messages passed. Includes the schedule index j + ncol(a.schedule$forward)
        message.names.mat <- rbind(message.names.mat, c(j+back.pass.off, a.schedule$backward[i, j], a.schedule$backward[i, j+1], msg.symb))
      }

    }
  }

  # A little reformatting for message.names.mat
  message.names.mat <- data.frame(as.numeric(message.names.mat[,1]), message.names.mat[,c(2:4)])
  colnames(message.names.mat) <- c("pass.num", "start.node", "end.node", "msg.symb")

  # At this point some messages are repeated because they came in on multiple paths to
  # the root. When doing sequential message passing, a node can't send out a message
  # until it has received all its incoming messages. So a messages pass number in the schedule
  # should correspond to the latest time it shows up in all the pass sequences. We can get that by
  # choosing the max pass number for a message in message.names.mat. Thats what we do now:
  unique.messages          <- unique(message.names.mat[,4])
  message.names.mat.pruned <- NULL

  for(i in 1:length(unique.messages)) {
    message.idxs          <- which(message.names.mat[,4] == unique.messages[i])
    message.names.mat.sub <- message.names.mat[message.idxs,]
    #print(message.names.mat.sub)

    message.pass.nums <- message.names.mat.sub[,1]
    #print(message.pass.nums)

    # Choose the last appearance of the message,
    # as it should have received all incoming messages
    # at that point. If there are many of the same max, just grab the first.
    max.message.pass.idx <- which(message.pass.nums == max(message.pass.nums))[1]
    #print(max.message.pass.idx)

    message.names.mat.pruned <- rbind(message.names.mat.pruned, message.names.mat.sub[max.message.pass.idx,])

  }
  message.names.mat <- message.names.mat.pruned
  message.names.mat <- message.names.mat[order(message.names.mat[,1]), ] # Re-order rows by scheduled pass number

  #print(message.names.mat)

  # Re-number the pass numbers to be consecutive in case they are not:
  pass.nums <- sort(unique(message.names.mat[,1]))
  for(i in 1:length(pass.nums)){
    message.names.mat[which(message.names.mat[,1] == pass.nums[i]), 1] <- i
  }
  #print(message.names.mat)

  # Container to hold all the messages passed over the network
  a.mailroom        <- rep(list(NULL), nrow(message.names.mat))
  names(a.mailroom) <- message.names.mat[,4]

  rownames(message.names.mat)   <- NULL
  message.container.info        <- list(a.mailroom, message.names.mat)
  names(message.container.info) <- c("message.container", "message.schedule.mat")

  return(message.container.info)

}


#' Get message type depending on name
#'
#' Get message type depending on name
#' a.name can me the message name or the name of the first node in the message
#'
#' The function will XXXX
#'
#' @param a.name a.name can me the message name or the name of the first node in the messageThe XX
#' @return The function will XX
#'
#' @export
message.type <- function(a.name) {

  # a.name can me the message name or the name of the first node in the message
  f.nodeQ <- ("f" ==  unlist(strsplit(x = a.name, split = ""))[1])

  if(f.nodeQ == T){
    mtyp <- "f2v"
  } else {
    mtyp <- "v2f"
  }

  return(mtyp)

}


#' Get all neighbors of start.node except end.node
#'
#' Get all neighbors of start.node except end.node
#' XXXX
#'
#' The function will XXXX
#'
#' @param a.name a.name can me the message name or the name of the first node in the messageThe XX
#' @return The function will XX
#'
#' @export
nex <- function(a.graphnel.obj, start.node, end.node) {

  # The neighbors:
  t.nes <- adj(a.graphnel.obj, start.node)[[1]]
  t.nes <- t.nes[-which(t.nes == end.node)]

  return(t.nes)

}



#' Grab incoming messages from the neighboring nodes and put into a list
#'
#' Grab incoming messages from the neighboring nodes and put into a list
#' Needed for constructing outgoing messages. NOTE: neighbor node set are all neighbors of
#' the start.node EXCEPT the end.node.
#'
#' The function will XXXX
#'
#' @param a.name a.name can me the message name or the name of the first node in the messageThe XX
#' @return The function will XX
#'
#' @export
get.incoming.messages <- function(start.node, end.node, factorgraph, message.list) {

  neibs             <- nex(factorgraph, start.node, end.node) # neighbors of start node except end node
  incoming.messages <- paste0(neibs, ".", start.node)         # incoming messages to start node
  #print(incoming.messages)

  message.names         <- names(message.list)
  message.idxs          <- sapply(1:length(incoming.messages), function(xx){which(message.names == incoming.messages[xx])})
  #print(message.idxs)

  if(length(message.idxs[[1]]) == 0){   # The starting node was a leaf, so no incoming messages
    # No names for this list indicates that the start node had no neighbors (except the end node)
    incoming.message.list <- list(NULL)
  } else {
    incoming.message.list <- message.list[message.idxs]
  }

  return(incoming.message.list)

}
