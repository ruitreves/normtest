#' Performs Wilcoxon/Mann-Whitney test on each row of a data.frame
#' @param my_data a numeric data.frame
#' @param var1 a list which maps columns of my_data to experimental factors
#' @param padj boolean, whether or not to include Benjamini Hochberg correction. defaults to FALSE
#' @return a data.frame of p values
#' @export
run_wilcox <- function(my_data, var1, padj = FALSE) {
    wilcox_res <- NULL
    for (i in 1:nrow(my_data)) {
        x <- suppressWarnings(stats::wilcox.test(unlist(my_data[i, ]) ~ var1, exact = TRUE))
        wilcox_res[[i]] <- x$p.value

        if (i%%100==0) {
            per <- round((i / nrow(my_data)) * 100)
            cat(paste("\r", "Wilcox-test ", per, "% done", sep=""))
        }
        else if (i == nrow(my_data)) {
            cat(paste("\r", "Wilcox-test ", "100% done", sep = ""))
        }
    }
    wilcox_res <- do.call(rbind, wilcox_res)
    wilcox_res <- as.data.frame(wilcox_res)
    rownames(wilcox_res) <- rownames(my_data)
    colnames(wilcox_res) <- "Pval"
    if (padj == TRUE) {
        wilcox_padj <- as.data.frame(stats::p.adjust(wilcox_res[, 1], method = "BH"))
        wilcox_res <- cbind(wilcox_res, wilcox_padj)
        colnames(wilcox_res) <- c("Pval", "Padj")
    }
    wilcox_res$test <- "wilcox-test"
    cat("\n")

    return(wilcox_res)
}
