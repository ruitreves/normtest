#' The "main" statomatic function. Takes an StatomaticDataset (sds) and performs all applicable statomatic analysis 
#' @param object An object of class "StatomaticDataset"
#' @param condition1 a list which maps columns of get_counts(sds) to experimental factors. If not specified, defaults to NULL and will use the design 
#' formula of sds
#' @param condition2 same as condition1. May be specified for a multi-factor experiment design, but defaults to NULL and if design formula is only of 
#' length 1, this condition is not used. 
#' @param padj boolean, whether or not to include Benjamini Hochberg correction. defaults to FALSE
#' @return a description of the method
#' @export 
setGeneric("sds_analyze", function(object, condition1 = NULL, condition2 = NULL, padj = FALSE) standardGeneric("sds_analyze"))
#'@inheritParams sds_analyze
#' @return a StatomaticDataset
#' @export
setMethod("sds_analyze", "StatomaticDataset", function(object, condition1 = NULL, condition2 = NULL, padj = FALSE) {
    if (is.null(SummarizedExperiment::assays(object)$counts)) stop("No counts found to analyze")
    counts <- SummarizedExperiment::assays(object)$counts
    f2 <- NULL
    if (is.null(condition1) && is.null(condition2)) {
        if (object@design == ~1) stop("You must define at least one condition to use for the model")
        condition1 <- all.vars(object@design)[1]
        f1 <- SummarizedExperiment::colData(object)[[condition1]]
        if (length(all.vars(object@design)) > 1) {
            condition2 <- all.vars(object@design)[2]
            f2 <- SummarizedExperiment::colData(object)[[condition2]]
        }
    }
    else if (is.null(condition1) && !is.null(condition2)) {
        if (condition2 %in% colnames(colData(object))) {
            f1 <- SummarizedExperiment::colData(object)[[condition2]]
        }
        else (stop(paste0(condition2, " not found in colData")))
    }
    else if (!is.null(condition1) && is.null(condition2)) {
        if (condition1 %in% colnames(colData(object))) {
            f1 <- SummarizedExperiment::colData(object)[[condition1]]
        }
        else (stop(paste0(condition1, " not found in colData")))
    }
    else {
        if (condition1 %in% colnames(colData(object)) && condition2 %in% colnames(colData(object))) {
            f1 <- SummarizedExperiment::colData(object)[[condition1]]
            f2 <- SummarizedExperiment::colData(object)[[condition2]]
        }
        else {(stop("The specified conditions were not found in colData"))}
    }

    var_length <- length(unique(f1))

    if (is.null(f2) && var_length < 3) {
        res <- twogroup_main(x = counts, var1 = f1, padj = padj, write_files = FALSE, include_test_name = TRUE)

        ne <- res[["ne"]]
        nu <- res[["nu"]]
        nn <- res[["nn"]]
        t_res <- res[[1]]
        welch <- res[[2]]
        wilcox <- res[[3]]
        fcs <- res[[4]]

        all_results <- rbind(t_res, welch, wilcox)

        object@results$t_test_results <- t_res
        object@results$welch_test_results <- welch
        object@results$wilcox_test_results <- wilcox
    }
    else {

        if (!is.null(f2)) {
            f3 <- SummarizedExperiment::colData(object)[, ncol(SummarizedExperiment::colData(sds))]
            res <- multigroup_main(x = counts, var1 = f1, var2 = f2, var3 = f3, padj = padj, write_files = FALSE)
            ne <- res[["ne"]]
            nu <- res[["nu"]]
            nn <- res[["nn"]]
            anova_results <- res[[1]]
            welch_results <- res[[2]]
            kw_results <- res[[3]]
            fcs <- res[[4]]
            if (padj == TRUE) {
                colnames(anova_results)[1:6] <- c(substitute(condition1), substitute(condition2), "interaction", 
                                        paste0(substitute(condition1), "_padj"), paste0(substitute(condition2), "_padj"), "interaction_padj")
                tukey_res <- anova_results[, -c(1:6)]
                dunnett_res <- welch_results[, -c(1:2)]
                dunn_res <- kw_results[, -c(1:2)]
            }
            else {
                colnames(anova_results)[1:3] <- c(substitute(condition1), substitute(condition2), "interaction")
                tukey_res <- anova_results[, -c(1:3)]
                dunnett_res <- welch_results[, -1]
                dunn_res <- kw_results[, -1]
            }
            all_results <- rbind(tukey_res, dunnett_res, dunn_res)
        }

        else {
            res <- multigroup_main(x = counts, var1 = f1, padj = padj, write_files = FALSE)
            ne <- res[["ne"]]
            nu <- res[["nu"]]
            nn <- res[["nn"]]
            anova_results <- res[[1]]
            welch_results <- res[[2]]
            kw_results <- res[[3]]
            fcs <- res[[4]]
            if (padj == TRUE) {
                colnames(anova_results)[1:2] <- c("Pval", "Padj")
                tukey_res <- anova_results[, -c(1:2)]
                dunnett_res <- welch_results[, -c(1:2)]
                dunn_res <- kw_results[, -c(1:2)]
            }
            else {
                colnames(anova_results)[1] <- "Pval"
                tukey_res <- anova_results[, -1]
                dunnett_res <- welch_results[, -1]
                dunn_res <- kw_results[, -1]
            }
            all_results <- rbind(anova_results, welch_results, kw_results)
        }

        object@results$anova_results <- anova_results
        object@results$welch_test_results <- welch_results
        object@results$kw_results <- kw_results
        object@results$tukey_test_results
        object@results$dunnett_test_results
        object@results$dunn_test_results
    }

    object@results$all_results <- all_results
    object@results$ne <- ne
    object@results$nu <- nu
    object@results$nn <- nn
    object@results$foldchange <- fcs

    return(object)
})