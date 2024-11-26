d <- dummy_data(group_max = 2)
x <- d[[1]]
sample_info <- d[[2]]

t_res <- run_ttest(x, sample_info$groups)

test_that("t_test works", {
  for (i in 1:nrow(x)) {
    test <- t.test(unlist(x[i, ]) ~ sample_info$groups, var.equal = TRUE)
    expect_equal(test$p.value, t_res[i, 1])

  }
})