

## construct_CovBlk #####

expect_equal(
  construct_CovBlk(c(3,3,3),c(0,0,0)),
  diag(3) *9
)
expect_equal(
  construct_CovBlk(c(3,3,3),c(1,2,3)),
  matrix(c(10,2,3,2,13,6,3,6,18),nrow=3)
)
expect_equal(
  construct_CovBlk(c(3,3,3),c(1,1,2),c(1,2,3)),
  matrix(c(11,3,5,3,14,8,5,8,22),nrow=3)
)

## construct_CovMat #####



# expect_equal_to_reference()
