#' Performs Kruskal-Wallis test on each row of a data.frame and reports the p value
#' @param my_data a numeric data.frame
#' @param var1 a list which maps columns of my_data to experimental factors
#' @param padj boolean, whether or not to include Benjamini Hochberg correction. defaults to FALSE
#' @return a data.frame of p values
#' @export
run_kruskal <- function(my_data, var1, padj = FALSE) {
    
    krusk_res <- NULL
    for (i in 1:nrow(my_data)) {
        #assign null each iteration
        krusk <- NULL
        #run test 
        krusk <- stats::kruskal.test(unlist(my_data[i, ]) ~ var1, as.data.frame(my_data))
        krusk <- krusk$p.value

        krusk_res <- rbind(krusk_res, krusk)

        if (i%%50==0) {
            per <- round((i / nrow(my_data)) * 100)
            cat(paste("\r", "Kruskal-Wallis Test ", per, "% done...", sep=""))
        }
        else if (i == nrow(my_data)) {
            cat(paste("\r", "Kruskal-Wallis Test ", "100% done", sep = ""))
        }
    }

    rownames(krusk_res) <- rownames(my_data)
    colnames(krusk_res) <- "Pval"

    if (padj == TRUE) {
        k_padj <- as.data.frame(stats::p.adjust(krusk_res[, 1], method = "BH"))
        krusk_res <- cbind(krusk_res, k_padj)
        colnames(krusk_res) <- c("Pval", "Padj")
    }

    cat("\n")

    return(krusk_res)
}