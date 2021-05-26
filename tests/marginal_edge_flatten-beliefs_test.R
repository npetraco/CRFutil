
dimnames(AB.bels[["Bel(A|B)"]])
dimnames(AB.bels[["Bel(B|A)"]])

dimnames(AB.bels[["Bel(A,B)"]])
AB.bels$`Bel(A,B)`
as.numeric(AB.bels$`Bel(A,B)`)
"col row"
"B=1 A=1"
"B=1 A=2"
"B=2 A=1"
"B=2 A=2"

dimnames(AB.bels[["Bel(A|B)"]])
AB.bels$`Bel(A|B)`
as.numeric(AB.bels$`Bel(A|B)`)
"col row"
"A=1 B=1"
"A=1 B=2"
"A=2 B=1"
"A=2 B=2"


dimnames(AB.bels[["Bel(B|A)"]])
AB.bels[["Bel(B|A)"]]
as.numeric(AB.bels$`Bel(B|A)`)
"col row"
"B=1 A=1"
"B=1 A=2"
"B=2 A=1"
"B=2 A=2"

flatten.marginal.edge.beliefs(AB.bels)

mij  <- AB.bels[["Bel(A,B)"]]
dnij <- dimnames(AB.bels[["Bel(A,B)"]])
#dnij
nij <- names(dnij)
#nij
#mij
for(i in 1:2){  # col
  for(j in 1:2){ # row
    #print(paste0("Col: ", nij[2], "=", i, " ", "Row: ",nij[1],"=",j, " Val: ", mij[j,i]))
    print(paste0(nij[2], "=", i, ",",nij[1],"=",j))
  }
}
AB.bels[["Bel(A|B)"]]

flatten.marginal.edge.beliefs(AB.bels)
AB.bels[["Bel(B)"]]


mij  <- AB.bels[["Bel(A)"]]
dnij <- dimnames(AB.bels[["Bel(A)"]])
#dnij
nij <- names(dnij)
#nij
#mij
for(i in 1:2){  # col
    #print(paste0("Col: ", nij[2], "=", i, " ", "Row: ",nij[1],"=",j, " Val: ", mij[j,i]))
    print(paste0(nij[1], "=", i))
}
