

## construct_trtMat #####
expect_equal(
  construct_timeadjust(Cl=c(1,1,1), timepoints=2, "factor"),
  matrix(c(rep(1,6),rep(0:1,3)),ncol=2)
)


## construct_timeadjust #####
expect_equal(
  construct_trtMat(Cl=c(1,1,1), trtDelay=NULL, dsntype="SWD", timepoints=4),
  { a <- matrix(0,nrow=3,ncol=4)
    a[upper.tri(a)] <- 1
    a
  }
)
expect_equal(
  construct_trtMat(Cl=c(1,1,1), trtDelay=NULL, dsntype="SWD", timepoints=5),
  { a <- matrix(0,nrow=3,ncol=5)
  a[upper.tri(a)] <- 1
  a
  }
)


  ## construct_DesMat #####
