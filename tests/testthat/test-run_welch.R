#mutli group test

d_multi <- dummy_data(group_min = 3)
x_multi <- d_multi[[1]]
sample_info_multi <- d_multi[[2]]

welch_res_multi <- run_welch(x_multi, sample_info_multi$group)
test_that("multi group welch works", {
for (i in 1:nrow(x_multi)) {
  w <- onewaytests::welch.test(x_multi[i, ] ~ sample_info_multi$group, data = as.data.frame(x_multi[i, ]), verbose = FALSE)
  expect_equal(w$p.value, welch_res_multi[i, 1])
}
})

#two group test 

d_tg <- dummy_data(group_max = 2)
x_tg <- d_tg[[1]]
sample_info_tg <- d_tg[[2]]

welch_res_tg <- run_welch(x_tg, sample_info_tg$groups)
test_that("two group welch works", {
  for (i in 1:nrow(x_tg)) {
    w <- onewaytests::welch.test(x_tg[i, ] ~ sample_info_tg$groups, data = as.data.frame(x_tg[i, ]), verbose = FALSE)
    expect_equal(w$p.value, welch_res_tg[i, 1])
  }
})