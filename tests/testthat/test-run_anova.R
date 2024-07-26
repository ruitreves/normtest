d <- dummy_data()
x <- d[[1]]
sample_info <- d[[2]]
sample_info$other <- c(rep("other1", ceiling((nrow(sample_info)) / 2)), rep("other2", floor((nrow(sample_info) / 2))))

res <- run_anova(x, sample_info$groups, sample_info$other)

for (i in 1:nrow(x)) {
  mod <- lm(unlist(x[i, ]) ~ sample_info$groups * sample_info$other)
  anova_res <- anova(mod)
  test <- anova_res$Pr[1:3]
  test <- na.omit(test)
  test <- unclass(test)
  test_that("multiplication works", {
    expect_setequal(test, res[i, ])
  })
}