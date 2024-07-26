#' For at least 3 groups. Test and sort data by distribution and variance equality
#' @param my_data A numeric set of data such as list, double, integer, etc.
#' @param var1 A column from a sample_info table used to define which column names of my_data go with var1
#' @param var2 A column from a sample_info table used to define which column names of my_data go with var2
#' 
#' @return A list of three dataframes ne, nu, nn
#' @export

test_norm <- function(my_data, var1, var2) {
    #genes with normally distributed residuals and equal variances, suitable for anova
    ne <- data.frame()
    #genes with normally distributed residuals but unequal variances, suitable for welch test
    nu <- data.frame()
    #genes with non-normally distributed residuals. 
    nn <- data.frame()
    for (i in 1:nrow(my_data)) {

        name_list <- rownames(my_data)

        model <- stats::lm(unlist(my_data[i, ]) ~ var1 * var2, as.data.frame(my_data))
        resids <- model$residuals
        #we first test if the residuals have equal variances between groups
        
        shap <- stats::shapiro.test(resids)


        if (shap$p.value >= 0.05) {
            #if normal, we test for equal variances
            vari <- car::leveneTest(resids ~ var1 * var2, as.data.frame(resids))
            #if normal and equal variances
            if (vari$Pr[1] >= 0.05) {
                ne <- rbind(ne, name_list[i])
            }
            #if normal and unequal variances
            else if (vari$Pr[1] < 0.05) {
                nu <- rbind(nu, name_list[i])
            }
        }

        else if (shap$p.value < 0.05) {
            #if not normal, we dont test for equal variance 
            nn <- rbind(nn, name_list[i])
        }

        #this is a just progess tracker for the function
        if (i%%100==0) {
            per <- round((i / nrow(my_data)) * 100)
            cat(paste0("\r", "sorting data ", per, "% done..."))
        }
    }
    #extract geneids for indexing

    if (nrow(ne) > 0) {
        x <- ne[, 1]
        ne <- my_data[x, ]
    }
    if (nrow(nu) > 0) {
         y <- nu[, 1]
         nu <- my_data[y, ]
    }
    if (nrow(nn) > 0) {
        z <- nn[, 1]
        nn <- my_data[z, ]
    }
    res <- c()
    res[[1]] <- ne
    res[[2]] <- nu
    res[[3]] <- nn

    cat("\n")
    return(res)
}