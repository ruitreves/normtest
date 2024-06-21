#' This function puts together all the multigroup functions from the package in the correct order.
#' @param my_data a data object
#' @param var1 a list that points column names of my_data to the factors of the experiment
#' @param var2 a list that points column names of my_data to the factors of the experiment
#' @param var3 a list that points column names of my_data to the factors of the experiment
#' @export


rui <- function(my_data, var1, var2, var3) {

    test_norm(my_data, var1, var2)

    ########
    ########    tier1
    ########


    if (nrow(tier1) > 0) {
        my_anova <- run_anova(tier1, var1, var2)

        tukey_res <- run_tukey(tier1, var3)

        fc_tier1 <- fold_change(tier1, var3)

        tukey_res <- combine_fc_pvals(tukey_res, fc_tier1)
        colnames(my_anova) <- c(substitute(var1), substitute(var2), "interaction")
        anova_results <- cbind(my_anova, tukey_res)
        write.csv(anova_results, "anova_results.csv")
    }

    ########
    ########    tier2
    ########

    if (nrow(tier2 > 0)) {    
        check_means(tier2, var3)

        welch_res <- run_welch(tier2, var3)

        dunnett_res <- run_dunnett(tier2, var3)

        fc_tier2 <- fold_change(tier2, var3)

        dunnett_res <- combine_fc_pvals(dunnett_res, fc_tier2)

        welch_results <- cbind(welch_res, dunnett_res)

        write.csv(welch_results, "welch_results.csv")
    }

    ########
    ########    tier3
    ########

    

    if (nrow(tier3) > 0) {        
        krusk_res <- run_kruskal(tier3, var3)

        #dunn_res <- run_dunn(tier3, var3)

        #using altered dunn d.test()
        dunn_res <- run_dunn(tier3, var3)

        fc_tier3 <- fold_change(tier3, var3)

        dunn_res <- combine_fc_pvals(dunn_res, fc_tier3)

        KW_results <- cbind(krusk_res, dunn_res)

        write.csv(KW_results, "KW_results.csv")
    }

}