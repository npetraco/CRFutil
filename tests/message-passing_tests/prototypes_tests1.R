init.message.storage <- function(a.schedule){

  a.mailroom    <- NULL
  message.names <- NULL

  # Initialize forward messages:
  for(i in 1:nrow(a.schedule$forward)){
    for(j in 1:(ncol(a.schedule$forward) - 1)){

      if(is.na(a.schedule$forward[i, j+1])) {
        break()
      } else {
        #print(paste0("Message: ", a.schedule$forward[i, j], "-->", a.schedule$forward[i, j+1]))
        message.names <- c(message.names, paste0(a.schedule$forward[i, j], ".", a.schedule$forward[i, j+1]))
        a.mailroom    <- c(a.mailroom, list(NULL)) # NULL can be interperited as an id message if necessary
      }

    }
  }

  # Initialize backward messages:
  for(i in 1:nrow(a.schedule$backward)){
    for(j in 1:(ncol(a.schedule$backward) - 1)){

      if(is.na(a.schedule$backward[i, j+1])) {
        break()
      } else {
        #print(paste0("Message: ", a.schedule$backward[i, j], "-->", a.schedule$backward[i, j+1]))
        message.names <- c(message.names, paste0(a.schedule$backward[i, j], ".", a.schedule$backward[i, j+1]))
        a.mailroom    <- c(a.mailroom, list(NULL)) # NULL can be interperited as an id message if necessary
      }

    }
  }


  names(a.mailroom) <- message.names
  return(a.mailroom)

}

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

# Get all neighbors of start.node except end.node
nex <- function(a.graphnel.obj, start.node, end.node) {

  # The neighbors:
  t.nes <- adj(a.graphnel.obj, start.node)[[1]]
  t.nes <- t.nes[-which(t.nes == end.node)]

  return(t.nes)

}


# Collect together messages from the neighboring nodes as a list
collect.messages <- function(neighborx.nodes, start.node, message.list) {

  req.messages <- rep(length(neighborx.nodes), list(NULL))
  print(req.messages)

}
