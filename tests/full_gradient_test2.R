samp.siz <- 10000
etest <- numeric(length(fit$par.stat))
etest
infinf <- infer.exact(fit)

# Nodes E.hattheta.phi
for(n in 1:fit$n.nodes){
  for(s in 1:fit$n.states[n]){
    if(fit$node.par[n,s,1] > 0){
      etest[fit$node.par[n,s,1]] <- etest[fit$node.par[n,s,1]] + samp.siz*infinf$node.bel[n,s]
    }
  }
}

# Edges E.hattheta.phi
for(e in 1:nrow(fit$edges)){
  n1xx <- fit$edges[e,1]
  n2xx <- fit$edges[e,2]
  for(es1xx in 1:fit$n.states[n1xx]){
    for(es2xx in 1:fit$n.states[n2xx]){
      if(fit$edge.par[[e]][es1xx,es2xx,1] > 0){
        etest[fit$edge.par[[e]][es1xx,es2xx,1]] <- etest[fit$edge.par[[e]][es1xx,es2xx,1]]  + samp.siz*infinf$edge.bel[[e]][es1xx,es2xx]
      }
    }
  }
}



etest/samp.siz
fit$par.stat/samp.siz


infinf$node.bel
fit$node.par[,,1]
infinf$node.bel[1,]*samp.siz
infinf$node.bel[2,]*samp.siz
infinf$node.bel[3,]*samp.siz
etest

fit$par.stat
infinf$edge.bel[[1]][1,1]*samp.siz + infinf$edge.bel[[1]][2,2]*samp.siz

rowSums(infinf$node.bel * (fit$node.par>0)[,,])
sum(infinf$edge.bel[[1]] * (fit$edge.par[[1]]>0)[,,])
sum(infinf$edge.bel[[2]] * (fit$edge.par[[2]]>0)[,,])
sum(infinf$edge.bel[[3]] * (fit$edge.par[[3]]>0)[,,])
etest/samp.siz
fit$par.stat/samp.siz
