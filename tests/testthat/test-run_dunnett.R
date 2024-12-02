library(PMCMRplus)

set.seed(112224)

d <- dummy_data(group_max = 6)
x <- d[[1]]
sample_info <- d[[2]]
z <- run_dunnett(x, sample_info$group)
z <- z[, -ncol(z)]

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
  for (i in 1:nrow(res)) {
    for (j in 1:ncol(res)) {
      expect_equal(res[i, j], z[i, j], tolerance = 0.001)
    }
  }
})

#test AdunnettT3test

counts_multi <- matrix(
  sample(1:1, 1600, replace = TRUE), 
  nrow = 100, ncol = 16,
  dimnames = list(c(paste0("gene_", 1:100)),
                  c(paste0("sample_", 1:16)))
)

metadata <- data.frame(samples = colnames(counts_multi), condition = c(rep("group1", 4), rep("group2", 4), rep("group3", 4), rep("group4", 4)))

x <- run_dunnett(counts_multi, metadata$condition)
x <- x[, -ncol(x)]

test_that("altered test works", {
  for (i in 1:nrow(x)) {
    for (j in 1:ncol(x)) {
      expect_equal(1, x[i, j])
    }
  }
})