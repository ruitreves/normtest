#' For exactly two groups. Test and sort data by distribution and equality of variances
#' @param my_data a data object
#' @param var1 a list that points column names of my_data to the factors of the experiment
#' @export 


two_group_test_norm <- function(my_data, var1) {
    #genes with normally distributed residuals and equal variances, suitable for anova
    ne <- data.frame()
    #genes with normally distributed residuals but unequal variances, suitable for welch test
    nu <- data.frame()
    #genes with non-normally distributed residuals. 
    nn <- data.frame()
    for (i in 1:nrow(my_data)) {

        name_list <- rownames(my_data)

        model <- stats::lm(unlist(my_data[i, ]) ~ var1, as.data.frame(my_data), na.action = na.exclude)
        resids <- model$residuals

        var1 <- as.factor(var1)

        #we first test if the residuals are normally distributed
        shap <- stats::shapiro.test(resids)
        
        if (shap$p.value >= 0.05) {
            #if normal, we test for equal variance
            vari <- car::leveneTest(resids ~ var1)

            if (vari$Pr[1] >= 0.05) {
                #if normal and equal variance
                ne <- rbind(ne, name_list[i])
            }

            else if (vari$Pr[1] < 0.05) {
                #if normal and unequal variance
                nu <- rbind(nu, name_list[i])
            }
        }

        else if (shap$p.value < 0.05) {
            #if not normally distributed, we do not test for equal variance
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