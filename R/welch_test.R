#' Performs Welch t-test on each row of a data.frame
#' @param my_data a numeric data.frame
#' @param var1 a list which maps columns of my_data to experimental factors 
#' @param padj boolean, whether or not to include Benjamini Hochberg correction. defaults to FALSE
#' @return a data.frame of p values
#' @export
run_welch <- function (my_data, var1, padj = FALSE) {
    welch_res <- NULL
    for (i in 1:nrow(my_data)) {
        welch <- onewaytests::welch.test(unlist(my_data[i, ]) ~ 
            var1, as.data.frame(my_data[i, ]), verbose = F)
        welch_res[[i]] <- welch$p.value
        if (i%%100 == 0) {
            per <- round((i/nrow(my_data)) * 100)
            cat(paste("\r", "Welch test ", per, "% done...", 
                sep = ""))
        }
        else if (i == nrow(my_data)) {
            cat(paste("\r", "Welch test 100% done", sep = ""))
        }
    }
    welch_res <- do.call(rbind, welch_res)
    welch_res <- as.data.frame(welch_res)
    rownames(welch_res) <- rownames(my_data)
    colnames(welch_res) <- "Pval"
    if (padj == TRUE) {
        welch_padj <- as.data.frame(stats::p.adjust(welch_res[, 
            1], method = "BH"))
        welch_res <- cbind(welch_res, welch_padj)
        colnames(welch_res) <- c("Pval", "Padj")
    }
    welch_res$test <- "welch-test"
    cat("\n")
    return(welch_res)
}


