#' Sorts data based on the distribution of it's residuals, as calculated by glm(family = gaussian), using shapiro and levene tests
#' @param my_data a numeric data.frame
#' @param var1 a list mapping columns of my_data to experimental factors
#' @param var2 same as var1, defaults to NULL
#' @return a list containing three data.frames, which are subsets of my_data. Access with double brackets. 
#' [[1]] is data which has normally distributed residuals and equal variance between groups, [[2]] has normally distributed residuals
#' and unequal variance between groups, [[3]] has non-normally distributed residuals. 
#' @export
test_norm <- function(my_data, var1, var2 = NULL) {
    if (is.null(var1)) stop("You must define at least one condition to use for the model")
    ne <- data.frame()
    nu <- data.frame()
    nn <- data.frame()

    name_list <- rownames(my_data)

    if (is.null(var2)) {
        formu1 <- as.formula(paste("unlist(my_data[i, ]) ~", "var1"))
        formu2 <- as.formula(paste("resids ~", "var1"))
    }
    else {
        formu1 <- as.formula(paste("unlist(my_data[i, ]) ~", "var1", "*", "var2"))
        formu2 <- as.formula(paste("resids ~", "var1", "*", "var2"))
    }

    for (i in 1:nrow(my_data)) {

        my_data <- as.data.frame(my_data)
        var1 <- as.factor(var1)

        model <- stats::glm(formu1, family = gaussian, data = my_data, na.action = na.exclude)
        resids <- model$residuals

        shap <- stats::shapiro.test(resids)
        
        if (shap$p.value >= 0.05) {
            vari <- car::leveneTest(formu2)

            if (vari$Pr[1] >= 0.05) {
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

