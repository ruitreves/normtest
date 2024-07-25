

d <- dummy_data()
x <- d[[1]]
sample_info <- d[[2]]

tukey_res <- run_tukey(x, sample_info$groups)

mod <- aov(unlist(x[1, ]) ~ sample_info$groups)
tt <- TukeyHSD(mod)
res <- as.data.frame(tt[1])
res <- t(res)
res <- as.data.frame(res[4, ])
res <- t(res)

big_res <- NULL
for (i in 1:nrow(x)) {
    mod <- aov(unlist(x[i, ]) ~ sample_info$groups)
    tt <- TukeyHSD(mod)
    res <- as.data.frame(tt[1])
    res <- t(res)
    res <- as.data.frame(res[4, ])
    res <- t(res)
    big_res <- rbind(big_res, res)
}

test_that("tukey works", {
  for (i in 1:nrow(big_res)) {
    for (j in 1:ncol(big_res)) {
        expect_equal(big_res[i, j], tukey_res[i, j])
    }
  }
})