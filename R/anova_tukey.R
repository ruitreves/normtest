#' Runs a two way anova based on var1 * var2 for each row of my_data
#' @param my_data a data object
#' @param var1 a list that points column names of my_data to the factors of the experiment
#' @param var2 a list that points column names of my_data to the factors of the experiment
#' @param padj logical, whether to add p adjust values or not. Calculated by p.adjust function with method = "BH"
#' @return a data object containing p values for var1, var2, and interaction of var1 * var2 for each row of my_data, computed by anova
#' @export

run_anova <- function(my_data, var1, var2, padj = FALSE) {
    my_anova <- NULL
    for (i in 1:nrow(my_data)) {
        #res <- t(as.data.frame(stats::anova(stats::lm(unlist(my_data[i,]) ~ var1 * var2))[1:3,5]))
        mod <- stats::lm(unlist(my_data[i, ]) ~ var1 * var2)
        a <- stats::anova(mod)
        res <- t(as.data.frame(a[1:3, 5]))
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
    colnames(my_anova) <- c(substitute(var1), substitute(var2), "interaction")
    rownames(my_anova) <- rownames(my_data)

    if (padj == TRUE) {
        var1_padj <- as.data.frame(stats::p.adjust(my_anova[, 1], method = "BH"))
        var2_padj <- as.data.frame(stats::p.adjust(my_anova[, 2], method = "BH"))
        intx_padj <- as.data.frame(stats::p.adjust(my_anova[, 3], method = "BH"))
        my_anova <- cbind(my_anova, var1_padj, var2_padj, intx_padj)
        colnames(my_anova) <- c(substitute(var1), substitute(var2), "interaction",
                                paste0(deparse(substitute(var1)), "_padj"), paste0(deparse(substitute(var2)), "_padj"),
                                "interaction_padj")
        my_anova <- my_anova[, c(1, 4, 2, 5, 3, 6)]
    }

    cat("\n")
    return(my_anova)
}

#'Runs a Tukey multicomparison test on each row of my_data based on var1
#' @param my_data a data object
#' @param var1 a list that points column names of my_data to the factors of the experiment
#' return a data object containing p values for each groupwise comparison as specified by var1, computed by TukeyHSD
#' @export

run_tukey <- function(my_data, var1) {
    test_tukey <- NULL
    for (i in 1:nrow(my_data)) {
        #my_data[i, ] <- unlist(my_data[i, ])
        model <- stats::aov(unlist(my_data[i, ]) ~ var1)
        test_result <- stats::TukeyHSD(model, conf.level=.95)
        res <- as.data.frame(test_result[1])
        res <- t(res)
        res <- as.data.frame(res[4,])
        res <- t(res) 
        test_tukey <- rbind(test_tukey, res)
        rownames(test_tukey)[i] <- rownames(my_data)[i]

        if (i%%100==0) {
            per <- round((i / nrow(my_data)) * 100)
            cat(paste("\r", "Tukey test ", per, "% done...", sep=""))
        }
        else if (i == nrow(my_data)) {
            cat(paste("\r", "Tukey test ", "100% done", sep = ""))
        }
    }
    x <- colnames(test_tukey)
    s <- strsplit(x, split = "-")
    srev <- c()
    for (i in 1:length(s)) {
        srev[[i]] <- rev(s[[i]])
        srev[[i]] <- paste(srev[[i]], collapse = "_vs_")
    }
    y <- paste("Pval_", srev, sep = "")
    colnames(test_tukey) <- y

    cat("\n")
    return(test_tukey)
}