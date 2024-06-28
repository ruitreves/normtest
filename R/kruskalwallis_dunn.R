#'Run Kruskal-Wallis test on each row of my_data
#' @param my_data a data object
#' @param var1 a list that points column names of my_data to the factors of the experiment
#' @return a data object with a p value for each row of my_data, computed using kruskal.test
#' @export

#Kruskal-Wallis test
run_kruskal <- function(my_data, var1) {
    
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
    colnames(krusk_res) <- "KW_Pval"

    cat("\n")

    return(krusk_res)
}

#'Runs Dunn's multi comparison test on each row of my_data
#' @param my_data a data object
#' @param var1 a list that points column names of my_data to the factors of the experiment
#' @return a data object with a pvalue for each groupwise comparison, as specified by var1, computed using an altered version of dunn.test
#' @export

#Dunn multi comparison test
run_dunn <- function(my_data, var1) {
    dunn_res <- NULL
    for (i in 1:nrow(my_data)) {
        #run dunn test
        dunn <- Adunn.test(unlist(my_data[i, ]), g=var1, method="bonferroni", table = FALSE, kw = FALSE, altp = TRUE)
        #we only want the 4th element, the padj value
        dunn_padj <- as.data.frame(dunn[4])
        #transpose 
        dunn_padj <- t(dunn_padj)
        #name columns as their comparisons
        #colnames(dunn) <- comp
        #add to final output
        dunn_res <- rbind(dunn_res, dunn_padj)

        if (i%%50==0) {
            per <- round((i / nrow(my_data)) * 100)
            cat(paste("\r", "Dunn Test ", per, "% done...", sep=""))
        }
        else if (i == nrow(my_data)) {
            cat(paste("\r", "Dunn Test ", "100% done", sep = ""))
        }
    }

    grp_list <- unique(var1)
    grp_list <- grp_list[order(grp_list)]
    col_name_list <- c()
    counter <- 1
    for (grp in grp_list) {
        counter <- counter + 1
        if (counter <= length(grp_list)) {
            for (j in counter:length(grp_list)) {
                the_name <- paste0(grp, " - ", grp_list[j])
                col_name_list <- append(col_name_list, the_name)
            }
        }
    }

    rownames(dunn_res) <- rownames(my_data)
    
    dunn_res <- dunn_res[, match(col_name_list, dunn$comparisons)]

    colnames(dunn_res) <- col_name_list
    colnames(dunn_res) <- gsub(" - ", "_vs_", colnames(dunn_res))

    cat("\n")

    return(dunn_res)
}