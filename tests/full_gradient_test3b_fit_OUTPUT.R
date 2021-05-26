# > w.test
# [1] -0.4689827 -0.2557522  0.1457067  0.8164156 -0.5966361  0.7967794  0.8893505  0.3215956  0.2582281 -0.8764275
# [11] -0.5880509 -0.6468865

#------------------------------------------------------------------------------------------


# Before negloglik(w.test, fit, samps, infer.exact)

# node/edge pots all 1

# > fit$par.stat <- mrf.stat(fit, samps)
# > fit$par.stat
# [1] 42 11 10 41  6  3 39  6 10 39  2  1
# > fit$gradient
# [1] 0 0 0 0 0 0 0 0 0 0 0 0
# > fit$nll
# [1] 0
# > fit$par
# [1] 0 0 0 0 0 0 0 0 0 0 0 0

#------------------------------------------------------------------------------------------

# After negloglik(w.test, fit, samps, infer.exact)    # NOTE, negloglik run here WAS MISSING - SIGN in from of w!!!!!*********
# > negloglik(w.test, fit, samps, infer.exact)
# [1] "num samps: 50"
# [1] "logZ: 2.77258872223978"
# [1] "offset: 138.629436111989"
# [1] "param vec:"
#  [1] -0.4689827 -0.2557522  0.1457067  0.8164156 -0.5966361  0.7967794  0.8893505  0.3215956  0.2582281 -0.8764275
# [11] -0.5880509 -0.6468865
# [1] "suff stats"
# [1] 42 11 10 41  6  3 39  6 10 39  2  1
# [,1]
# [1,] 153.0524

# node/edge pots STILL all 1 *********

# > fit$par.stat
# [1] 42 11 10 41  6  3 39  6 10 39  2  1
# > fit$gradient
# [1] 0 0 0 0 0 0 0 0 0 0 0 0
# > fit$nll
# [1] 0
# > fit$par
# [1] 0 0 0 0 0 0 0 0 0 0 0 0

# When infer.exact is run immediately afterward on  fit object we get: *******

# > infer.exact(fit)
# $node.bel
# [,1] [,2]
# [1,]  0.5  0.5
# [2,]  0.5  0.5
# [3,]  0.5  0.5
# [4,]  0.5  0.5
#
# $edge.bel
# $edge.bel[[1]]
# [,1] [,2]
# [1,] 0.25 0.25
# [2,] 0.25 0.25
#
# $edge.bel[[2]]
# [,1] [,2]
# [1,] 0.25 0.25
# [2,] 0.25 0.25
#
# $edge.bel[[3]]
# [,1] [,2]
# [1,] 0.25 0.25
# [2,] 0.25 0.25
#
# $edge.bel[[4]]
# [,1] [,2]
# [1,] 0.25 0.25
# [2,] 0.25 0.25
#
#
# $logZ
# [1] 2.772589



# Update potentials with w.test:
# > upd.pots[[1]]
# , , 1
#
# [,1] [,2]
# [1,] 0.6256384    1
# [2,] 0.7743338    1
# [3,] 1.1568569    1
# [4,] 2.2623760    1
#
# > upd.pots[[2]]
# [[1]]
# [,1]     [,2]
# [1,] 0.5506609 1.000000
# [2,] 1.0000000 2.218385
#
# [[2]]
# [,1]     [,2]
# [1,] 2.433549 1.000000
# [2,] 1.000000 1.379327
#
# [[3]]
# [,1]      [,2]
# [1,] 1.294634 1.0000000
# [2,] 1.000000 0.4162674
#
# [[4]]
# [,1]      [,2]
# [1,] 0.5554088 1.0000000
# [2,] 1.0000000 0.5236737


# Above pots Substituted in fit object and run infer exact:
# > infer.exact(fit)
# $node.bel
# [,1]      [,2]
# [1,] 0.3419521 0.6580479
# [2,] 0.3799990 0.6200010
# [3,] 0.6355479 0.3644521
# [4,] 0.6567970 0.3432030
#
# $edge.bel
# $edge.bel[[1]]
# [,1]      [,2]
# [1,] 0.1427452 0.1992069
# [2,] 0.2372537 0.4207942
#
# $edge.bel[[2]]
# [,1]       [,2]
# [1,] 0.2801022 0.06184992
# [2,] 0.3766948 0.28135310
#
# $edge.bel[[3]]
# [,1]      [,2]
# [1,] 0.2067635 0.1732355
# [2,] 0.4287844 0.1912166
#
# $edge.bel[[4]]
# [,1]       [,2]
# [1,] 0.3587479 0.27680005
# [2,] 0.2980491 0.06640297
#
#
# $logZ
# [1] 3.119087

# When Re-run negloglik with the fit object now (pots from w.test have been substituted):
# > negloglik(w.test, fit, samps, infer.exact)
# [1] "num samps: 50"
# [1] "logZ: 3.11908661275298"
# [1] "offset: 155.954330637649"
# [1] "param vec:"
# [1] -0.4689827 -0.2557522  0.1457067  0.8164156 -0.5966361  0.7967794  0.8893505  0.3215956  0.2582281 -0.8764275
# [11] -0.5880509 -0.6468865
# [1] "suff stats"
# [1] 42 11 10 41  6  3 39  6 10 39  2  1
# [,1]
# [1,] 141.5314


# When negloglik is run MAKE SURE:
# 1. pots in fit obj correspond to w.test sent in.
# They will probably need to be updated. Needed to get correct Z from infer function.
# Think this was causing much of the problem.
# ****NOTE:
#  negative log likelihood IS  -w %*% suffStat + nInstances * infer.info$logZ
#
# 2. par in fit object is updated with w.test ??????
#
# 3. nll in fit object is updated
