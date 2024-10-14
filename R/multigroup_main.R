#' This function puts together all the multigroup functions from the package in the correct order.
#' @param my_data a data object
#' @param var1 a list that points column names of my_data to the factors of the experiment
#' @param var2 a list that points column names of my_data to the factors of the experiment
#' @param var3 a list that points column names of my_data to the factors of the experiment
#' @export


multigroup_main <- function(my_data, var1, var2, var3, padj = FALSE) {

    res <- test_norm(my_data, var1, var2)
    
    ne <- res[[1]]
    nu <- res[[2]]
    nn <- res[[3]]

    ########
    ########    ne
    ########
    

    if (nrow(ne) > 0) {
        my_anova <- run_anova(ne, var1 = var1, var2 = var2, padj)

        if (padj == TRUE) {
            colnames(my_anova) <- c(substitute(var1), paste0(deparse(substitute(var1)), "_padj"), substitute(var2),
                                    paste0(deparse(substitute(var2)), "_padj"), "interaction", "interaction_padj")
        }
        else {
            colnames(my_anova) <- c(substitute(var1), substitute(var2), "interaction")
        }

        tukey_res <- run_tukey(ne, var3)

        fc_ne <- fold_change(ne, var3)

        tukey_res <- combine_fc_pvals(tukey_res, fc_ne)
        anova_results <- cbind(my_anova, tukey_res)
        write.csv(anova_results, "anova_results.csv")
    }

    ########
    ########    nu
    ########

    if (nrow(nu > 0)) {    
        check_means(nu, var3)

        welch_res <- run_welch(nu, var3, padj)

        dunnett_res <- run_dunnett(nu, var3)

        fc_nu <- fold_change(nu, var3)

        dunnett_res <- combine_fc_pvals(dunnett_res, fc_nu)

        welch_results <- cbind(welch_res, dunnett_res)

        write.csv(welch_results, "multi_group_welch_results.csv")
    }

    ########
    ########    nn
    ########

    if (nrow(nn) > 0) {        
        krusk_res <- run_kruskal(nn, var3, padj)

        #dunn_res <- run_dunn(nn, var3)

        #using altered dunn d.test()
        dunn_res <- run_dunn(nn, var3)

        fc_nn <- fold_change(nn, var3)

        dunn_res <- combine_fc_pvals(dunn_res, fc_nn)

        KW_results <- cbind(krusk_res, dunn_res)

        write.csv(KW_results, "KW_results.csv")
    }
}