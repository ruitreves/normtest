#' Performs a t-test on each row of a data.frame. 
#' @param my_data a numeric data.frame
#' @param var1 a list which maps columns of my_data to experimental factors. This function can handle exactly two experimental factos
#' @param padj boolean, whether or not to include Benjamini Hochberg correction. defaults to FALSE
#' @return a data.frame of p values
#' @export
run_ttest <- function (my_data, var1, padj = FALSE) {
    t_res <- NULL
    for (i in 1:nrow(my_data)) {
        x <- stats::t.test(unlist(my_data[i, ]) ~ var1, var.equal = TRUE)
        t_res[[i]] <- x$p.value
        if (i%%100 == 0) {
            per <- round((i/nrow(my_data)) * 100)
            cat(paste("\r", "t-test ", per, "% done", sep = ""))
        }
        else if (i == nrow(my_data)) {
            cat(paste("\r", "t-test ", "100% done", sep = ""))
        }
    }
    t_res <- do.call(rbind, t_res)
    t_res <- as.data.frame(t_res)
    rownames(t_res) <- rownames(my_data)
    colnames(t_res) <- "Pval"
    if (padj == TRUE) {
        t_padj <- as.data.frame(stats::p.adjust(t_res[, 1], method = "BH"))
        t_res <- cbind(t_res, t_padj)
        colnames(t_res) <- c("Pval", "Padj")
    }
    t_res$test <- "t-test"
    cat("\n")
    return(t_res)
}


