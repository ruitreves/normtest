#' This function puts together all the two group functions from the package in the correct order.
#' @param my_data a data object
#' @param var1 a list that points column names of my_data to the factors of the experiment
#' @export

two_group_rui <- function(my_data, var1) {
    #sort data into normal and equal variance, normal and nonequal variances,
    #and nonnormal.
    two_group_test_norm(my_data, var1)

    if (nrow(tier1) > 0) {
        t_res <- run_ttest(tier1, var1)
        tier1_folds <- fold_change(tier1, unique(var1))

        t_res <- cbind(t_res, tier1_folds)
        write.csv(t_res, "t_test_results.csv")
    }

    if (nrow(tier2) > 0) {

        welch_res <- run_welch(tier2, var1)
        tier2_folds <- fold_change(tier2, unique(var1))

        welch_res <- cbind(welch_res, tier2_folds)
        write.csv(welch_res, "welch_results.csv")
    }

    if (nrow(tier3) > 0) {
        wilcox_res <- run_wilcox(tier3, var1)
        tier3_folds <- fold_change(tier3, unique(var1))

        wilcox_res <- cbind(wilcox_res, tier3_folds)
        write.csv(wilcox_res, "wilcox_results.csv")
    }
    cat("Done.")
}
