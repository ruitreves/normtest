#' The fold_change functoin iterates through the groups specified by var1
#' to calculate the log2FoldChange of mean(groupA) / mean(groupB)
#' @param my_data a dataframe, double, list with numeric entries
#' @param var1 a list that points column names of my_data to the factors of the experiment
#' @return a dataframe containing foldchanges for each groupwise comparison
#' @export

fold_change <- function(my_data, var1) {
    cat("Calculating fold changes...")
    big_res <- NULL
    grp_list <- unique(var1)
    grp_list <- grp_list[order(grp_list)]
    for (i in rownames(my_data)) {
        x <- my_data[i, ]
        small_res <- NULL
        counter <- 1

        for (grp in grp_list) {
            counter <- counter + 1
            mean1 <- mean(unlist(x)[var1 == grp])

            if (counter <= length(grp_list)) {
                for (j in counter:length(grp_list)) {
                    mean2 <- mean(unlist(x)[var1 == grp_list[j]])

                    if (mean2 == 0) {
                        fld_cng <- Inf
                    }
                    else if (mean1 == 0) {
                        fld_cng <- -Inf
                    }
                    else {
                        av <- abs(mean1 / mean2)
                        fld_cng <- log2(av)
                    }

                    small_res <- cbind(small_res, fld_cng)
                }
            }
        }
        big_res <- rbind(big_res, small_res)
    }

    col_name_list <- c()
    counter <- 1
    for (grp in grp_list) {
        counter <- counter + 1
        if (counter <= length(grp_list)) {
            for (j in counter:length(grp_list)) {
                the_name <- paste("log2FoldChange_", grp_list[j], "_vs_", grp, sep="")
                col_name_list <- append(col_name_list, the_name)
            }
        }
    }

    rownames(big_res) <- rownames(my_data)
    colnames(big_res) <- col_name_list
    cat("\n")
    return(big_res)
}



#' combine fold changes and pvals into single dataframe. takes column
#' i from data1 and column i from data2 and puts them next to each other
#' @param data1 a data object with n columns and m rows
#' @param data2 a data object with n columns and m rows
#' @return a dataframe made of data1 and data2
#' @export
combine_fc_pvals <- function(data1, data2) {
    #placeholders
    out_df <- data.frame()
    list_names <- c()
    #dataframes should have the same dimensions
    for (i in 1:ncol(data1)) {
        #bind first columns of inputs, then second columns, etc... 
        x <- cbind(as.data.frame(data1[, i]), as.data.frame(data2[, i]))
        #append each colname to list_names
        list_names <- append(list_names, colnames(data1)[i])
        list_names <- append(list_names, colnames(data2)[i])
        #we rbind the transpose of x, which we will transpose again later
        out_df <- rbind(out_df, t(x))
    }
    #transpose again
    out_df <- t(out_df)
    #rename cols
    colnames(out_df) <- list_names

    return(out_df)
}

#' check if two or more group means are zero
#' @param my_data a data object 
#' @param var1 a list that points column names of my_data to the factors of the experiment
#' @export
check_means <- function(my_data, var1) {
    grp_list <- unique(var1)
    bad_list <- c()
    for (i in 1:nrow(my_data)) {    
        counter <- 1
        x <- my_data[i, ]

        for (grp in grp_list) {
            counter <- counter + 1
            mean1 <- mean(unlist(x)[var1 == grp])

            if (counter <= length(grp_list)) {
                for (j in counter:length(grp_list)) {
                    mean2 <- mean(unlist(x)[var1 == grp_list[j]])

                    if (mean1 == 0 & mean2 ==0) {
                        bad_list <- append(bad_list, i)
                    }
                }
            }
        }
    }

    if (length(bad_list) > 0) {
        my_data <- my_data[-bad_list, ]
        print(paste(length(bad_list), "rows removed from", substitute(my_data)))
        assign("tier2", my_data, envir = parent.frame())
    }
}