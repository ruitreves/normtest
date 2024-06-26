#'Runs a Welch t-test on each row of my_data
#' @param my_data a data object
#' @param var1 a list that points column names of my_data to the factors of the experiment
#' @return a data object with a p value for each row of my_data, computed using welch.test
#' @export

#Welch test
run_welch <- function(my_data, var1) {
    welch_res <- NULL
    for (i in 1:nrow(my_data)) {
        welch <- onewaytests::welch.test(unlist(my_data[i, ]) ~ var1, as.data.frame(my_data[i,]), verbose=F)
        welch_pval <- welch$p.value
        welch_res <- rbind(welch_res, welch_pval)

        #for style points, 
        if (i%%100==0) {
            per <- round((i / nrow(my_data)) * 100)
            cat(paste("\r", "Welch test ", per, "% done...", sep=""))
        }
        else if (i == nrow(my_data)) {
            cat(paste("\r", "Welch test 100% done", sep = ""))
        }
    }

    rownames(welch_res) <- rownames(my_data)
    colnames(welch_res) <- "Pval"
    welch_res <- as.data.frame(welch_res)

    cat("\n")

    return(welch_res)
}

#'Runs DunnettT3 multi comparison test on each row of my_data
#' @param my_data a data object
#' @param var1 a list that points column names of my_data to the factors of the experiment
#' @return a data object containing p values for each groupwise comparison as specified by var1 
#' @export
#DunnettT3  multi comparison test

run_dunnett <- function(my_data, var1) {

    dunnett_res <- NULL
    for (i in 1:nrow(my_data)) {
        dunnet <- PMCMRplus::dunnettT3Test(unlist(my_data[i, ]), g = as.factor(var1))
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
                colnames(pvals)[counter] <- paste(rn[b], cn[a], sep = "-vs-")
                counter <- counter + 1
            }
        }
        pvals <- pvals[, colSums(is.na(pvals)) == 0]
        dunnett_res <- rbind(dunnett_res, pvals)

        if (i%%70==0) {
            per <- round((i / nrow(my_data)) * 100)
            cat(paste("\r", "Dunnett T3 test ", per, "% done...", sep=""))
        }
        else if (i == nrow(my_data)) {
            cat(paste("\r", "Dunnett T3 test 100% done", sep = ""))
        }
    }
    rownames(dunnett_res) <- rownames(my_data)
    x <- colnames(dunnett_res)
    x <- paste("Pval_", x, sep = "")
    colnames(dunnett_res) <- x

    cat("\n")

    return(dunnett_res)
}