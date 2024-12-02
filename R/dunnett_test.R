#' Performs Dunnett's T3 test on each row of a data.frame
#' @param my_data a numeric data.frame
#' @param var1 a list which maps columns of my_data to experimental factors
#' @return a data.frame of p values
#' @export

run_dunnett <- function(my_data, var1) {

    dunnett_res <- NULL
    for (i in 1:nrow(my_data)) {
        dunnet <- AdunnettT3Test(unlist(my_data[i, ]), g = as.factor(var1))
        dunnet <- dunnet$p.value

        rn <- rownames(dunnet)
        cn <- colnames(dunnet)

        pvals <- NULL
        for (n in dunnet) {
            pvals <- cbind(pvals, n)
        }

        counter <- 1
        for (a in 1:length(cn)) {
            for (b in 1:length(rn)) {
                colnames(pvals)[counter] <- paste(cn[a], rn[b], sep = "_vs_")
                counter <- counter + 1
            }
        }
        pvals <- pvals[, colSums(is.na(pvals)) == 0]
        dunnett_res <- rbind(dunnett_res, pvals)
        dunnett_res <- as.data.frame(dunnett_res)

        if (i%%70==0) {
            per <- round((i / nrow(my_data)) * 100)
            cat(paste("\r", "Dunnett T3 test ", per, "% done...", sep=""))
        }
        else if (i == nrow(my_data)) {
            cat(paste("\r", "DunnettT3 test 100% done", sep = ""))
        }
    }
    rownames(dunnett_res) <- rownames(my_data)
    colnames(dunnett_res) <- paste0("Pval_", colnames(dunnett_res))
    dunnett_res$test <- "dunnett-test"
    cat("\n")

    return(dunnett_res)
}

