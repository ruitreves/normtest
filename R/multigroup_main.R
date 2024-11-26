#' The multi-group statomatic procedure, pulling together all applicable functions in the correct order. Use this function to analyze a dataset 
#' with more than 2 groups
#' @param x a numeric data.frame
#' @param var1 a list which maps columns of x to experimental factors. Must be exactly 2 factors
#' @param var2 same as var1, defaults to NULL
#' @param var3 same as var2
#' @param padj boolean, whether or not to include Benjamini Hochberg correction. defaults to FALSE
#' @param write_files boolean, whether or not to write csv files of the test results instead of returning. Defaults to FALSE
#' @param file_names a list of three file names to be used when writing files. Specify as a vector using c(). The first will be the name of the t-test results,
#' the second the welch-test results, and the third the wilcox-test results. Only used if write_files = TRUE
#' @return a list containing all test results, fold changes, and subsets of x as sorted by test_norm, unless write_files = TRUE, in which case nothing is returned.
#' @export
multigroup_main <- function(x, var1, var2 = NULL, var3 = NULL, padj = FALSE, write_files = FALSE, file_names = c("anova_results.csv", "welch_results.csv", "kw_results.csv")) {
    if (write_files == TRUE && length(file_names) != 3) stop("This function creates three files but has not been provided three names")
    if (is.null(var3)) {var3 <- var1}
    nenunn <- test_norm(x, var1, var2)
    ne <- nenunn[[1]]
    nu <- nenunn[[2]]
    nn <- nenunn[[3]]
    all_fcs <- NULL
    anova_results <- data.frame()
    welch_results <- data.frame()
    kw_results <- data.frame()
    if (nrow(ne) > 0) {
        my_anova <- run_anova(ne, var1 = var1, var2 = var2, padj)
        if (padj == TRUE) {
            if (is.null(var2)) {
                colnames(my_anova) <- c(substitute(var1), paste0(deparse(substitute(var1)), "_padj"))
            }
            else {
                colnames(my_anova) <- c(substitute(var1), substitute(var2), "interaction", 
                                    paste0(deparse(substitute(var1)), "_padj"), paste0(deparse(substitute(var2)), "_padj"), "interaction_padj")
            }
        }
        else {
            if (is.null(var2)) {
                colnames(my_anova) <- "Pval"
            }
            else {
                colnames(my_anova) <- c(substitute(var1), substitute(var2), "interaction")
            }
            
        }
        tukey_res <- run_tukey(ne, var3)
        fc_ne <- fold_change(ne, var3)
        tukey_res <- combine_fc_pvals(tukey_res, fc_ne)
        anova_results <- cbind(my_anova, tukey_res)
        all_fcs <- rbind(all_fcs, fc_ne)
    }
    if (nrow(nu) > 0) {
        check_means(nu, var3)
        welch_res <- run_welch(nu, var3, padj)
        dunnett_res <- run_dunnett(nu, var3)
        fc_nu <- fold_change(nu, var3)
        dunnett_res <- combine_fc_pvals(dunnett_res, fc_nu)
        welch_results <- cbind(welch_res, dunnett_res)
        all_fcs <- rbind(all_fcs, fc_nu)
    }
    if (nrow(nn) > 0) {
        krusk_res <- run_kruskal(nn, var3, padj)
        dunn_res <- run_dunn(nn, var3)
        fc_nn <- fold_change(nn, var3)
        dunn_res <- combine_fc_pvals(dunn_res, fc_nn)
        kw_results <- cbind(krusk_res, dunn_res)
        all_fcs <- rbind(all_fcs, fc_nn)
    }

    if (write_files == TRUE) {
        write.csv(anova_results, file_names[1])
        write.csv(welch_results, file_names[2])
        write.csv(kw_results, file_names[3])
        return(NULL)
    }
    
    res <- list("anova_results" = anova_results, "welch_results" = welch_results, "kw_results" = kw_results, "foldchanges" = all_fcs,
                "ne" = ne, "nu" = nu, "nn" = nn)
    return(res)
}