test_that("welch works", {
library(onewaytests)
new_dat <- data.frame()
cols <- c()
for (i in 1:20) {
    n <- rnorm(1000, 0, 1)
    x <- paste("X", i, sep="")
    y <- paste("Y", i, sep="")
    cols <- append(cols, y)
    new_dat <- rbind(new_dat, n)
}

new_dat <- t(new_dat)

rows <- c()
for (i in 1:1000) {
    x <- paste("X", i, sep="")
    rows <- append(rows, x)
}

rownames(new_dat) <- rows
colnames(new_dat) <- cols

sample_info <- data.frame(sample = colnames(new_dat), group = c(rep("Group1", 5), rep("Group2", 5), rep("Group3", 5), rep("Group4", 5)))

welch_res <- run_welch(new_dat, sample_info$group)

for (i in 1:nrow(new_dat)) {
  w <- welch.test(new_dat[i, ] ~ sample_info$group, data = as.data.frame(new_dat[i, ]), verbose = FALSE)
  expect_equal(w$p.value, welch_res[i, ])
}
})

