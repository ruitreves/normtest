#' Performs Tukey's HSD test on each row of a data.frame
#' @param my_data a numeric data.frame
#' @param var1 a list which maps columns of my_data to experimental factors
#' @return a data frame of p values
#' @export
run_tukey <- function(my_data, var1) {
    res <- c()
    for (i in 1:nrow(my_data)) {
        model <- stats::aov(unlist(my_data[i, ]) ~ var1)
        test_result <- stats::TukeyHSD(model, conf.level=.95)
        test_result <- as.data.frame(test_result[1])
        test_result <- t(test_result[, 4, drop = FALSE])
        res[[i]] <- test_result

        if (i%%100==0) {
            per <- round((i / nrow(my_data)) * 100)
            cat(paste0("\r", "Tukey test ", per, "% done..."))
        }
        else if (i == nrow(my_data)) {
            cat(paste0("\r", "Tukey test ", "100% done"))
        }
    }
    test_tukey <- do.call(rbind, res)
    test_tukey <- as.data.frame(test_tukey)
    x <- colnames(test_tukey)
    s <- strsplit(x, split = "-")
    srev <- c()
    for (i in 1:length(s)) {
        srev[[i]] <- rev(s[[i]])
        srev[[i]] <- paste(srev[[i]], collapse = "_vs_")
    }
    y <- paste("Pval_", srev, sep = "")
    colnames(test_tukey) <- y
    rownames(test_tukey) <- rownames(my_data)
    test_tukey$test <- "tukey-test"
    cat("\n")
    return(test_tukey)
}


