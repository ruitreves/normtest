d <- dummy_data(group_max = 2)
x <- d[[1]]
sample_info <- d[[2]]

res <- run_wilcox(x, sample_info$groups)

test_that("wilcox works", {
  for (i in 1:nrow(x)) {
    w <- wilcox.test(unlist(x[i, ]) ~ sample_info$groups, exact = TRUE)
    expect_equal(w$p.value, res[i, ])
  }
})
