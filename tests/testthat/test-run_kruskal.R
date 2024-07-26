
d <- dummy_data()
x <- d[[1]]
sample_info <- d[[2]]

krusk_res <- run_kruskal(x, sample_info$groups)

for (i in 1:nrow(x)) {
  krusk <- stats::kruskal.test(unlist(x[i, ]) ~ sample_info$groups, as.data.frame(x))
  krusk <- krusk$p.value
  test_that("multiplication works", {
    expect_equal(krusk, krusk_res[i, ])
  })
}