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
        assign("nu", my_data, envir = parent.frame())
    }
}

#' Set rownames to a specified column and remove that column from the dataframe
#' @param dat a data.frame
#' @param column_index the index of the column you want as your rownames, defaults to 1
#' @return the original data frame with the row names set as specified and the specified column removed
#' @export
cf <- function(dat, column_index = 1) {
    rownames(dat) <- dat[, column_index]
    dat[, column_index] <- NULL
    return(dat)
}

#' Take n elements of a list and concatenate them. Will introduce NAs if length(x) %% n != 0
#' @param x a list
#' @param n an integer specifying the number of elements to concatenate together
#' @param sep a string, the separating string between elements once concatenated
#' @export

takeN <- function(x, n = 1, sep = "") {
    counter <- 1
    for (i in 1:(ceiling(length(x) / n))) {
        for (j in counter:(counter + (n-1))) {
            cat(x[j])
            cat(sep)
        }
        cat("\n")
        counter <- counter + n        
    }
}

#' Chop multi group result files into single comparisons.
#' @param x a data frame, usually one produced by running the main multigroup normtest function, containing only groupwise comparisons.
#' @param p p value cutoff for significance. Defaults to 0.05
#' @param f log2 fold change cut off. Defaults to 1
#' @export 

chop <- function(x, p = 0.05, f = 1) {
    chopped <- c()
    res <- c()
    counter <- 1
    for (i in 1:(ncol(x) / 2)) {
        pfc <- x[, counter:(counter + 1)]
        res[[i]] <- pfc
        counter <- counter + 2
    }
    sig_res <- c()
    for (j in 1:length(res)) {
        cur <- res[[j]]
        sigs <- cur[cur[, 1] < p, ]
        sigs <- sigs[abs(sigs[, 2]) > f, ]
        sig_res[[j]] <- sigs
    }
    #chopped <- c(res, sig_res)
    chopped[[1]] <- res
    chopped[[2]] <- sig_res
    return(chopped)
}

#' Produces a dummy data set from a normal distribution and sample_info object. 
#' @param group_min Minimum random number of groups. Defaults to 3
#' @param group_max Maximum random number of groups. Defaults to 20. Note: group_max * dim max must be <= 1000
#' @param dim_min Minimum random dimension of data frame. Defaults to 1
#' @param dim_max Maximum random dimension of data frame. Defaults to 50. Note: group_max * dim max must be <= 1000
#' @param param_min Minimum random parameters for normal distribution (i.e., mean and variance). Defaults to 0
#' @param param_max Maximum random parameters for normal distribution. Defaults to 20
#' @return A list, element one is the randomly generated data frame, element two is sample_info. Access with double brackets: [[1]], [[2]].
#' @export 

dummy_data <- function(group_min = 2, group_max = 20, dim_min = 1, dim_max = 50, param_min = 0, param_max = 20) {
    group_num <- floor(stats::runif(1, min = group_min, max = group_max))

    group_names <- c()
    for (i in 1:group_num) {
        group_name <- paste0("group", i)
        group_names <- append(group_names, group_name)
    }

    dims <- ceiling(stats::runif(2, min = dim_min, max = dim_max))
    col_num <- dims[1] * group_num
    row_num <- dims[2]

    parms <- ceiling(stats::runif(2, min = param_min, max = param_max))
    m <- parms[2]
    v <- parms[1]

    x <- data.frame()
    cols <- c()
    for (i in 1:col_num) {
        n <- rnorm(row_num, m, v)
        c <- paste0("Y", i)
        cols <- append(cols, c)
        x <- rbind(x, n)
    }

    x <- t(x)

    rows <- c()
    for (i in 1:row_num) {
        r <- paste0("X", i)
        rows <- append(rows, r)
    }

    rownames(x) <- rows
    colnames(x) <- cols

    sample_info <- data.frame(samples = colnames(x), groups = group_names)

    dum <- c()
    dum[[1]] <- x
    dum[[2]] <- sample_info

    return(dum)
}
