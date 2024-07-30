d <- dummy_data(group_min = 3)
x <- d[[1]]
sample_info <- d[[2]]

tukey_res <- run_tukey(x, sample_info$groups)

test_that("tukey 2", {
  for (i in 1:nrow(x)) {
    mod <- aov(unlist(x[i, ]) ~ sample_info$groups)
    test <- TukeyHSD(mod)
    test1 <- test[1]
    test_res <- as.data.frame(test1)
    for (j in 1:length(test_res[, 4])) {
      expect_equal(test_res[, 4][j], tukey_res[i, ][[j]])
    }
  }
})
