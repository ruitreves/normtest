d <- dummy_data()
x <- d[[1]]
sample_info <- d[[2]]

welch_res <- run_welch(x, sample_info$group)
test_that("welch works", {
for (i in 1:nrow(x)) {
  w <- onewaytests::welch.test(x[i, ] ~ sample_info$group, data = as.data.frame(x[i, ]), verbose = FALSE)
  expect_equal(w$p.value, welch_res[i, ])
}
})

