#' Runs an unpaired t-test on each row of my_data and reports the p value, computed by t.test
#' @param my_data a data object
#' @param var1 list that points column names of my_data to the factors of the experiment
#' @param padj logical, whether to add p adjust values or not. Calculated by p.adjust function with method = "BH"
#' @export 

run_ttest <- function(my_data, var1, padj = FALSE) {
    t_res <- NULL
    for(i in 1:nrow(my_data)) {
        x <- stats::t.test(unlist(my_data[i, ]) ~ var1, var.equal = TRUE)
        t_res <- rbind(t_res, x$p.value)

        if (i%%100==0) {
        per <- round((i / nrow(my_data)) * 100)
        cat(paste("\r", "t-test ", per, "% done", sep=""))
        }
        else if (i == nrow(my_data)) {
            cat(paste("\r", "t-test ", "100% done", sep = ""))
        }
    }

    rownames(t_res) <- rownames(my_data)
    colnames(t_res) <- "Pval"

    if (padj == TRUE) {
        t_padj <- as.data.frame(stats::p.adjust(t_res[, 1], method = "BH"))
        t_res <- cbind(t_res, t_padj)
        colnames(t_res) <- c("Pval", "Padj")
    }

    cat("\n")

    return(t_res)
}
