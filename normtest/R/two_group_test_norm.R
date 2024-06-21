#' For exactly two groups. Test and sort data by distribution and equality of variances
#' @param my_data a data object
#' @param var1 a list that points column names of my_data to the factors of the experiment
#' @export 


two_group_test_norm <- function(my_data, var1) {
    #genes with normally distributed residuals and equal variances, suitable for anova
    tier1 <- data.frame()
    #genes with normally distributed residuals but unequal variances, suitable for welch test
    tier2 <- data.frame()
    #genes with non-normally distributed residuals. 
    tier3 <- data.frame()
    for (i in 1:nrow(my_data)) {

        name_list <- rownames(my_data)

        model <- stats::lm(unlist(my_data[i, ]) ~ var1, as.data.frame(my_data))
        resids <- model$residuals

        var1 <- as.factor(var1)

        #we first test if the residuals are normally distributed
        shap <- stats::shapiro.test(resids)
        
        if (shap$p.value >= 0.05) {
            #if normal, we test for equal variance
            vari <- car::leveneTest(resids ~ var1)

            if (vari$Pr[1] >= 0.05) {
                #if normal and equal variance
                tier1 <- rbind(tier1, name_list[i])
            }

            else if (vari$Pr[1] < 0.05) {
                #if normal and unequal variance
                tier2 <- rbind(tier2, name_list[i])
            }
        }

        else if (shap$p.value < 0.05) {
            #if not normally distributed, we do not test for equal variance
            tier3 <- rbind(tier3, name_list[i])
        }

        #this is a just progess tracker for the function
        if (i%%100==0) {
            per <- round((i / nrow(my_data)) * 100)
            cat(paste("\r", "sorting data ", per, "% done...", sep=""))
        }
    }
    #extract geneids for indexing
    x <- tier1[,1]
    y <- tier2[,1]
    z <- tier3[,1]

    tier1 <- my_data[x,]
    tier2 <- my_data[y,]
    tier3 <- my_data[z,]

    assign("tier1", tier1, envir = .GlobalEnv)
    assign("tier2", tier2, envir = .GlobalEnv)
    assign("tier3", tier3, envir = .GlobalEnv)
    cat("\n")
}