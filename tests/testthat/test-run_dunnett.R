library(PMCMRplus)

d <- dummy_data(group_max = 6)
x <- d[[1]]
sample_info <- d[[2]]
z <- run_dunnett(x, sample_info$group)

res <- NULL
for (i in 1:nrow(x)) {
  d <- PMCMRplus::dunnettT3Test(unlist(x[i, ]), g = as.factor(sample_info$group))
  p <- d$p.value

  pvals <- NULL
  for (i in p) {
      pvals <- cbind(pvals, i)
  }
  pvals <- pvals[, colSums(is.na(pvals)) == 0]
  pvals <- t(as.data.frame(pvals))

  res <- rbind(res, pvals)
}

colnames(res) <- colnames(z)
rownames(res) <- rownames(z)

test_that("dunnett works", {
  #skip("skip")
  for (i in 1:nrow(res)) {
    for (j in 1:ncol(res)) {
      expect_equal(res[i, j], z[i, j], tolerance = 0.001)
    }
  }
})
