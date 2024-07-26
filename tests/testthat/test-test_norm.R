d <- dummy_data()
x <- d[[1]]
sample_info <- d[[2]]

sample_info$other <- c(rep("other1", ceiling((nrow(sample_info)) / 2)), rep("other2", floor((nrow(sample_info) / 2))))

res <- test_norm(x, sample_info$groups, sample_info$other)

pass_shap <- NULL
for (i in 1:nrow(x)) {
    rows <- rownames(x)
    mod <- glm(unlist(x[i, ]) ~ sample_info$groups * sample_info$other)
    resids <- mod$residuals
    shap <- shapiro.test(resids)
    if (shap$p.value >= 0.05) {
        pass_shap <- rbind(pass_shap, rows[i])
    }
}
pass1 <- x[pass_shap, ]

library(car)
pass_lev <- NULL
for (i in 1:nrow(pass1)) {
    rows <- rownames(pass1)
    mod <- glm(unlist(pass1[i, ]) ~ sample_info$groups * sample_info$other)
    resids <- mod$residuals
    vari <- suppressWarnings(leveneTest(resids ~ sample_info$groups * sample_info$other))
    if (vari$Pr[1] >= 0.05) {
        pass_lev <- rbind(pass_lev, rows[i])
    }
}
pass2 <- x[pass_lev, ]

test_that("shaprio works", {
  expect_setequal(pass2, res[[1]])
})
test_that("levene works", {
  expect_setequal(setdiff(pass1, pass2), res[[2]])
})
test_that("both work", {
  expect_setequal(setdiff(x, union(pass1, pass2)), res[[3]])
})