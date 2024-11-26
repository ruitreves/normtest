#' Performs one or two way anova on each row of a data.frame and reports the p value
#' @param my_data a numeric data.frame
#' @param var1 a list which maps columns of my_data to experimental factors
#' @param var2 same as var1, but defaults to NULL. using var1 and var2 result in two way anova (if both are valid)
#' @param padj boolean, whether or not to include Benjamini Hochberg correction. defaults to FALSE
#' @return a data.frame of p values
#' @export
run_anova <- function(my_data, var1, var2 = NULL, padj = FALSE) {
    my_anova <- NULL
    if (is.null(var1)) stop("You must define at least one condition to use for the model")
    if (is.null(var2)) {
        formu <- as.formula("unlist(my_data[i, ]) ~ var1")
        idx <- 1
        cols <- "Pval"
        
        cat("Using one way anova \n")
    }
    else {
        formu <- as.formula("unlist(my_data[i, ]) ~ var1 * var2")
        idx <- 1:3
        cols <- c(substitute(var1), substitute(var2), "interaction")
        cat("Using two way anova \n")
    }

    for (i in 1:nrow(my_data)) {
        mod <- stats::lm(formu)
        a <- stats::anova(mod)
        res <- t(as.data.frame(a[idx, 5]))
        my_anova <- as.data.frame(rbind(my_anova, res))

        #this is a just progess tracker for the function
        if (i%%100==0) {
            per <- round((i / nrow(my_data)) * 100)
            cat(paste("\r", "anova ", per, "% done", sep=""))
        }
        else if (i == nrow(my_data)) {
            cat(paste("\r", "anova ", "100% done", sep = ""))
        }
    }

    colnames(my_anova) <- cols
    rownames(my_anova) <- rownames(my_data)
    if (padj == TRUE) {
        padj <- as.data.frame(lapply(my_anova, p.adjust, method = "BH"))
        colnames(padj) <- paste0("padj_", colnames(padj))
        my_anova <- cbind(my_anova, padj)
    }

    cat("\n")
    return(my_anova)
}


