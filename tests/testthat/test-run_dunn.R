d <- dummy_data()
x <- d[[1]]
sample_info <- d[[2]]

dunn_res <- run_dunn(x, sample_info$groups)

dunn_test <- NULL
for (i in 1:nrow(x)) {
    dunn <- Adunn.test(unlist(x[i, ]), g=sample_info$groups, method="bonferroni", table = FALSE, kw = FALSE, altp = TRUE)
    dunn_padj <- as.data.frame(dunn[4])
    dunn_padj <- t(dunn_padj)

    dunn_test <- rbind(dunn_test, dunn_padj)
}

grp_list <- unique(sample_info$groups)
grp_list <- grp_list[order(grp_list)]

col_names <- c()
counter <- 1
for (g in grp_list) {
    counter <- counter +  1
    if (counter <= length(grp_list)) {
        for (i in counter:length(grp_list)) {
        the_name <- paste0(g, " - ", grp_list[i])
        col_names <- append(col_names, the_name)
        }
    }
}

dunn_test <- dunn_test[, match(col_names, dunn$comparisons)]
colnames(dunn_test) <- col_names

test_that("dunn works", {
  for (i in 1:nrow(dunn_test)) {
    for (j in 1:ncol(dunn_test)) {
        expect_equal(dunn_test[i, j], dunn_res[i, j])
    }
  }
})
