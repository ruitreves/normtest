#'Runs two sided Wilcoxon Rank Sum test on each row of my_data and reports the p value, computed by wilcox.test
#' @param my_data a data object
#' @param var1 list that points column names of my_data to the factors of the experiment
#' @export 

run_wilcox <- function(my_data, var1) {
    wilcox_res <- NULL
    for (i in 1:nrow(my_data)) {
        x <- suppressWarnings(stats::wilcox.test(unlist(my_data[i, ]) ~ var1, exact = TRUE))
        wilcox_res <- rbind(wilcox_res, x$p.value)

        if (i%%100==0) {
            per <- round((i / nrow(my_data)) * 100)
            cat(paste("\r", "Wilcox-test ", per, "% done", sep=""))
        }
        else if (i == nrow(my_data)) {
            cat(paste("\r", "Wilcox-test ", "100% done", sep = ""))
        }
    }

    rownames(wilcox_res) <- rownames(my_data)
    colnames(wilcox_res) <- "Pval"

    cat("\n")

    return(wilcox_res)
}