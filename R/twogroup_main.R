#' This function puts together all the two group functions from the package in the correct order.
#' @param my_data a data object
#' @param var1 a list that points column names of my_data to the factors of the experiment
#' @export

twogroup_main <- function(my_data, var1, padj = FALSE) {
    #sort data into normal and equal variance, normal and nonequal variances,
    #and nonnormal.
    res <- two_group_test_norm(my_data, var1)
    ne <- res[[1]]
    nu <- res[[2]]
    nn <- res[[3]]

    if (nrow(ne) > 0) {
        t_res <- run_ttest(ne, var1, padj)
        ne_folds <- fold_change(ne, var1)

        t_res <- cbind(t_res, ne_folds)
        write.csv(t_res, "t_test_results.csv")
    }

    if (nrow(nu) > 0) {

        welch_res <- run_welch(nu, var1, padj)
        nu_folds <- fold_change(nu, var1)

        welch_res <- cbind(welch_res, nu_folds)
        write.csv(welch_res, "two_group_welch_results.csv")
    }

    if (nrow(nn) > 0) {
        wilcox_res <- run_wilcox(nn, var1, padj)
        nn_folds <- fold_change(nn, var1)

        wilcox_res <- cbind(wilcox_res, nn_folds)
        write.csv(wilcox_res, "wilcox_results.csv")
    }
    cat("Done.")
    cat("\n")
}
