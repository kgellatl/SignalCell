set.seed(0)
nrows <- 1000L
ncols <- 100L
vals <- sample(
  x = c(0, 1, 2),
  prob = c(0.8, 0.1, 0.1),
  size = nrows * ncols,
  replace = TRUE
)
matBaseR <- matrix(vals, nrow = nrows)
matBaseR
##      [,1] [,2] [,3] [,4] [,5] [,6]
## [1,]    2    1    0    0    0    0
## [2,]    0    0    0    0    0    1
## [3,]    0    2    0    0    1    0
## [4,]    0    1    0    0    0    0

# Convert to sparse format
matSparse <- as(matBaseR, "sparseMatrix")
matSparse
