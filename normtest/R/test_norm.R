#test and sort data by distribution
test_norm <- function(my_data, var1, var2) {
    #genes with normally distributed residuals and equal variances, suitable for anova
    tier1 <- data.frame()
    #genes with normally distributed residuals but unequal variances, suitable for welch test
    tier2 <- data.frame()
    #genes with non-normally distributed residuals. 
    tier3 <- data.frame()
    library(car)
    for (i in 1:nrow(my_data)) {

        name_list <- rownames(my_data)

        model <- lm(unlist(my_data[i, ]) ~ var1 * var2, as.data.frame(my_data))
        resids <- model$residuals
        #we first test if the residuals have equal variances between groups
        
        shap <- shapiro.test(resids)


        if (shap$p.value >= 0.05) {
            #if normal, we test for equal variances
            vari <- leveneTest(resids ~ var1 * var2, as.data.frame(resids))
            #if normal and equal variances
            if (vari$Pr[1] >= 0.05) {
                tier1 <- rbind(tier1, name_list[i])
            }
            #if normal and unequal variances
            else if (vari$Pr[1] < 0.05) {
                tier2 <- rbind(tier2, name_list[i])
            }
        }

        else if (shap$p.value < 0.05) {
            #if not normal, we dont test for equal variance 
            tier3 <- rbind(tier3, name_list[i])
        }

        #this is a just progess tracker for the function
        if (i%%100==0) {
            per <- round((i / nrow(my_data)) * 100)
            cat(paste("\r", "sorting data ", per, "% done...", sep=""))
        }
    }
    #extract geneids for indexing

    if (nrow(tier1) > 0) {
        x <- tier1[, 1]
        tier1 <- my_data[x, ]
    }
    if (nrow(tier2) > 0) {
         y <- tier2[, 1]
         tier2 <- my_data[y, ]
    }
    if (nrow(tier3) > 0) {
        z <- tier3[, 1]
        tier3 <- my_data[z, ]
    }
    assign("tier1", tier1, envir = .GlobalEnv)
    assign("tier2", tier2, envir = .GlobalEnv)
    assign("tier3", tier3, envir = .GlobalEnv)
    cat("\n")
}