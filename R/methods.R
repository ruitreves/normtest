#' Runs the test_norm function on get_counts(sds)
#' @param object A StatomaticDataset
#' @param condition1 a list which maps columns of get_counts(sds) to experimental factors. If not specified, defaults to NULL and will use the design 
#' formula of sds
#' @param condition2 same as condition1. May be specified for a multi-factor experiment design, but defaults to NULL and if design formula is only of 
#' length 1, this condition is not used. 
#' @return a description of the method 
setGeneric("sds_testnorm", function(object, condition1 = NULL, condition2 = NULL) standardGeneric("sds_testnorm"))
#' @inheritParams sds_testnorm
#' @return a StatomaticDataset
#' @export
setMethod("sds_testnorm", "StatomaticDataset", function(object, condition1 = NULL, condition2 = NULL) {
    if (is.null(assays(object)$counts)) stop("No counts found to analyze")
    f2 <- NULL
    if (is.null(condition1) && is.null(condition2)) {
        if (object@design == ~1) stop("You must define at least one condition to use for the model")
        design_var1 <- all.vars(object@design)[1]
        f1 <- colData(object)[[design_var1]]
        if (length(all.vars(object@design)) > 1) {
            design_var2 <- all.vars(object@design)[2]
            f2 <- colData(object)[[design_var2]]
        }
    }
    else if (is.null(condition1) && !is.null(condition2)) {
        if (condition2 %in% colnames(colData(object))) {
            f1 <- colData(object)[[condition2]]
        }
        else (stop(paste0(condition2, " not found in colData")))
    }
    else if (!is.null(condition1) && is.null(condition2)) {
        if (condition1 %in% colnames(colData(object))) {
            f1 <- colData(object)[[condition1]]
        }
        else (stop(paste0(condition1, " not found in colData")))
    }
    else {
        if (condition1 %in% colnames(colData(object)) && condition2 %in% colnames(colData(object))) {
            f1 <- colData(object)[[condition1]]
            f2 <- colData(object)[[condition2]]
        }
        else {(stop("The specified conditions were not found in colData"))}
    }
    res <- test_norm(my_data = assays(object)$counts, var1 = f1, var2 = f2)
    object@results$ne <- res[[1]]
    object@results$nu <- res[[2]]
    object@results$nn <- res[[3]]

    return(object)
})

#' Runs fold_change function on data contained in a StatomaticDataset
#' @param object A StatomaticDataset
#' @param condition The experimental factor to consider for calculating fold changes. Defaults to NULL and will use the design contained in sds
#' @return A description of the method
#'
setGeneric("sds_foldchange", function(object, condition = NULL) standardGeneric("sds_foldchange"))
#' @inheritParams sds_foldchange
#' @return a StatomaticDatset
#' @export
setMethod("sds_foldchange", "StatomaticDataset", function(object, condition = NULL) {
    if (is.null(assays(object)$counts)) stop("No counts found")
    if (is.null(condition)) {
        if (length(all.vars(object@design)) == 1) {
            design_var <- all.vars(object@design)[1]
            var1 <- colData(object)[[design_var]]
        }
        else if (length(all.vars(object@design)) > 1) {
            var1 <- colData(object)[[ncol(colData(object))]]
        }
        else {stop("Please specify a condition to compare")}
    }
    else {
        if (condition %in% colnames(colData(object))) {
            var1 <- colData(object)[[condition]]
        }
        else {stop(paste0("Condition: ", condition, " not found in colData"))}
    }

    object@results$foldchange <- fold_change(my_data = assays(object)$counts, var1 = var1)
    return(object)
})

#' Runs anova on data contained in a StatomaticDataset
#' @param object A StatomaticDataset
#' @param condition1 a list which maps columns of get_counts(sds) to experimental factors. If not specified, defaults to NULL and will use the design 
#' formula of sds
#' @param condition2 same as condition1. if specified here or in the sds design formula will result in a two way anova 
#' @param padj boolean, whether or not to include Benjamini Hochberg correction. defaults to FALSE
#' @return a description of the method
#' @export
setGeneric("sds_anova", function(object, condition1 = NULL, condition2 = NULL, padj = FALSE) standardGeneric("sds_anova"))
#' @inheritParams sds_anova
#' @return A StatomaticDataset
#' @export
setMethod("sds_anova", "StatomaticDataset", function(object, condition1 = NULL, condition2 = NULL, padj = FALSE) {
    if (is.null(assays(object)$counts)) stop("No counts found to analyze")
    f2 <- NULL
    if (is.null(condition1) && is.null(condition2)) {
        if (object@design == ~1) stop("You must define at least one condition to use for the model")
        design_var1 <- all.vars(object@design)[1]
        f1 <- colData(object)[[design_var1]]
        cols <- c(substitute(design_var1))
        f2 <- NULL
        if (length(all.vars(object@design)) > 1) {
            design_var2 <- all.vars(object@design)[2]
            f2 <- colData(object)[[design_var2]]
            cols <- c(substitute(design_var1), substitute(design_var2), "interaction")
        }
    }
    else if (is.null(condition1) && !is.null(condition2)) {
        if (condition2 %in% colnames(colData(object))) {
            f1 <- colData(object)[[condition2]]
            cols <- c(substitute(condition2))
        }
        else (stop(paste0(condition2, " not found in colData")))
    }
    else if (!is.null(condition1) && is.null(condition2)) {
        if (condition1 %in% colnames(colData(object))) {
            f1 <- colData(object)[[condition1]]
            cols <- c(substitute(condition1))
        }
        else (stop(paste0(condition1, " not found in colData")))
    }
    else {
        if (condition1 %in% colnames(colData(object)) && condition2 %in% colnames(colData(object))) {
            f1 <- colData(object)[[condition1]]
            f2 <- colData(object)[[condition2]]
            cols <- c(substitute(condition1), substitute(condition2), "interaction")
        }
        else {(stop("The specified conditions were not found in colData"))}
    }
    if (padj == TRUE) {
        cols <- append(cols, c(paste0("padj_", cols)))
    }
    res <- run_anova(assays(object)$counts, var1 = f1, var2 = f2, padj = padj)
    colnames(res) <- cols
    object@results$anova_results <- res
    return(object)
})

#' Runs tukey-test on data in a StatomaticDataset
#' @param object a StatomaticDataset
#' @param condition a list which maps columns of get_counts(sds) to experimental factors. If not specified, defaults to NULL and will use the design 
#' formula of sds
#' @return a description of the method
#' @export
setGeneric("sds_tukeytest", function(object, condition = NULL) standardGeneric("sds_tukeytest"))
#' @inheritParams sds_tukeytest
#' @return a StatomaticDatset
#' @export
setMethod("sds_tukeytest", "StatomaticDataset", function(object, condition = NULL) {
    if (is.null(assays(object)$counts)) stop("No counts found to analyze")

    if (is.null(condition)) {
        if (length(all.vars(object@design)) == 1) {
            design_var <- all.vars(object@design)[1]
            var1 <- colData(object)[[design_var]]
            cat("Tukey-test on", substitute(design_var), "\n")
        }
        else {
            design_var <- colnames(colData(object))[ncol(colData(object))]
            var1 <- colData(object)[[design_var]]
            cat("Tukey-test on", substitute(design_var), "\n")
        }
    }
    else {
        if (condition %in% colnames(colData(object))) {
            var1 <- colData(object)[[condition]]
            cat("Tukey-test on", substitute(condition), "\n")
        }
        else (stop(paste0("Condition: ", condition, " not found in colData")))
    }
    object@results$tukey_test_results <- run_tukey(my_data = assays(object)$counts, var1 = var1)
    return(object)
})

#' Run welch-test on data contained in a StatomaticDataset
#' @param object a StatomaticDataset
#' @param condition a list which maps columns of get_counts(sds) to experimental factors. If not specified, defaults to NULL and will use the design 
#' formula of sds
#' @param padj boolean, whether or not to include Benjamini Hochberg correction. defaults to FALSE
#' @return a description of the method
#' @export 
setGeneric("sds_welchtest", function(object, condition = NULL, padj = FALSE) standardGeneric("sds_welchtest"))
#' @inheritParams sds_welchtest
#' @return a StatomaticDataset
#' @export
setMethod("sds_welchtest", "StatomaticDataset", function(object, condition = NULL, padj = FALSE) {
    if (is.null(assays(object)$counts)) stop("No counts found to analyze")

    if (is.null(condition)) {
        if (length(all.vars(object@design)) == 1) {
            design_var <- all.vars(object@design)[1]
            var1 <- colData(object)[[design_var]]
            cat("Welch-test on", substitute(design_var), "\n")
        }
        else {
            design_var <- colnames(colData(object))[ncol(colData(object))]
            var1 <- colData(object)[[design_var]]
            cat("Welch-test on", substitute(design_var), "\n")
        }
    }
    else {
        if (condition %in% colnames(colData(object))) {
            var1 <- colData(object)[[condition]]
            cat("Welch-test on", substitute(condition), "\n")
        }
        else (stop(paste0("Condition: ", condition, " not found in colData")))
    }
    object@results$welch_test_results <- run_welch(my_data = assays(object)$counts, var1 = var1, padj = padj)
    return(object)
})

#' Runs Dunnett T3 test on data contained in a Statomatic Dataset
#' @param object a StatomaticDataset
#' @param condition a list which maps columns of get_counts(sds) to experimental factors. If not specified, defaults to NULL and will use the design 
#' formula of sds
#' @return a description of the method 
#' @export 
setGeneric("sds_dunnetttest", function(object, condition = NULL) standardGeneric("sds_dunnetttest"))
#' @inheritParams sds_dunnetttest
#' @return a StatomaticDataset
#' @export
setMethod("sds_dunnetttest", "StatomaticDataset", function(object, condition = NULL) {
    if (is.null(assays(object)$counts)) stop("No counts found to analyze")

    if (is.null(condition)) {
        if (length(all.vars(object@design)) == 1) {
            design_var <- all.vars(object@design)[1]
            var1 <- colData(object)[[design_var]]
            cat("Dunnett-test on", substitute(design_var), "\n")
        }
        else {
            design_var <- colnames(colData(object))[ncol(colData(object))]
            var1 <- colData(object)[[design_var]]
            cat("Dunnett-test on", substitute(design_var), "\n")
        }
    }
    else {
        if (condition %in% colnames(colData(object))) {
            var1 <- colData(object)[[condition]]
            cat("Dunnett-test on", substitute(condition), "\n")
        }
        else (stop(paste0("Condition: ", condition, " not found in colData")))
    }
    object@results$dunnett_test_results <- run_dunnett(my_data = assays(object)$counts, var1 = var1)
    return(object)
})

#' Runs Kruskal-Wallis test on data contained in a StatomaticDataset
#' @param object a StatomaticDataset
#' @param condition a list which maps columns of get_counts(sds) to experimental factors. If not specified, defaults to NULL and will use the design 
#' formula of sds
#' @param padj boolean, whether or not to include Benjamini Hochberg correction. defaults to FALSE
#' @return a description of the method
#' @export 
setGeneric("sds_kruskaltest", function(object, condition = NULL, padj = FALSE) standardGeneric("sds_kruskaltest"))
#' @inheritParams sds_kruskaltest
#' @return a StatomaticDataset
#' @export
setMethod("sds_kruskaltest", "StatomaticDataset", function(object, condition = NULL, padj = FALSE) {
    if (is.null(assays(object)$counts)) stop("No counts found to analyze")

    if (is.null(condition)) {
        if (length(all.vars(object@design)) == 1) {
            design_var <- all.vars(object@design)[1]
            var1 <- colData(object)[[design_var]]
            cat("Kruskal-wallis-test on", substitute(design_var), "\n")
        }
        else {
            design_var <- colnames(colData(object))[ncol(colData(object))]
            var1 <- colData(object)[[design_var]]
            cat("Kruskal-wallis-test on", substitute(design_var), "\n")
        }
    }
    else {
        if (condition %in% colnames(colData(object))) {
            var1 <- colData(object)[[condition]]
            cat("Kruskal-wallis-test on", substitute(condition), "\n")
        }
        else (stop(paste0("Condition: ", condition, " not found in colData")))
    }
    object@results$kw_test_results <- run_kruskal(my_data = assays(object)$counts, var1 = var1, padj = padj)
    return(object)
})

#' Runs Dunn test on data contained in a StatomaticDataset
#' @param object a StatomaticDataset
#' @param condition a list which maps columns of get_counts(sds) to experimental factors. If not specified, defaults to NULL and will use the design 
#' formula of sds
#' @return a description of the method
#' @export
setGeneric("sds_dunntest", function(object, condition = NULL) standardGeneric("sds_dunntest"))
#' @inheritParams sds_dunntest
#' @return a StatomaticDataset
#' @export
setMethod("sds_dunntest", "StatomaticDataset", function(object, condition = NULL) {
    if (is.null(assays(object)$counts)) stop("No counts found to analyze")

    if (is.null(condition)) {
        if (length(all.vars(object@design)) == 1) {
            design_var <- all.vars(object@design)[1]
            var1 <- colData(object)[[design_var]]
            cat("Dunn-test on", substitute(design_var), "\n")
        }
        else {
            design_var <- colnames(colData(object))[ncol(colData(object))]
            var1 <- colData(object)[[design_var]]
            cat("Dunn-test on", substitute(design_var), "\n")
        }
    }
    else {
        if (condition %in% colnames(colData(object))) {
            var1 <- colData(object)[[condition]]
            cat("Dunn-test on", substitute(condition), "\n")
        }
        else (stop(paste0("Condition: ", condition, " not found in colData")))
    }
    object@results$dunn_test_results <- run_dunn(my_data = assays(object)$count, var1 = var1)
    return(object)
})

#' Runs t-test on data contained in a StatomaticDataset
#' @param object a StatomaticDataset
#' @param padj boolean, whether or not to include Benjamini Hochberg correction. defaults to FALSE
#' @return a description of the method 
#' @export
setGeneric("sds_ttest", function(object, condition = NULL, padj = FALSE) standardGeneric("sds_ttest"))
#' @inheritParams sds_ttest
#' @return a StatomaticDataset
#' @export
setMethod("sds_ttest", "StatomaticDataset", function(object, condition = NULL, padj = FALSE) {
    if (is.null(assays(object)$counts)) stop("No counts found")

    if (is.null(condition)) {
        if (length(all.vars(object@design)) == 1) {
            design_var <- all.vars(object@design)[1]
            var1 <- colData(object)[[design_var]]
            cat("t-test on", substitute(design_var), "\n")
        }
        else {
            design_var <- colnames(colData(object))[ncol(colData(object))]
            var1 <- colData(object)[[design_var]]
            cat("t-test on", substitute(design_var), "\n")
        }
    }
    else {
        if (condition %in% colnames(colData(object))) {
            var1 <- colData(object)[[condition]]
            cat("t-test on", substitute(condition), "\n")
        }
        else (stop(paste0("Condition: ", condition, " not found in colData")))
    }

    object@results$t_test_results <- run_ttest(my_data = assays(object)$counts, var1 = var1, padj = padj)
    return(object)
})

#welch test compatible with multi and two group

#' Runs wilcox test on data contained in a StatomaticDataset
#' @param object a StatomaticDataset
#' @param condition a list which maps columns of get_counts(sds) to experimental factors. If not specified, defaults to NULL and will use the design 
#' formula of sds
#' @param padj boolean, whether or not to include Benjamini Hochberg correction. defaults to FALSE
#' @return a description of the method 
#' @export
setGeneric("sds_wilcoxtest", function(object, condition = NULL, padj = FALSE) standardGeneric("sds_wilcoxtest"))
#' @inheritParams sds_wilcoxtest
#' @return a StatomaticDataset
#' @export
setMethod("sds_wilcoxtest", "StatomaticDataset", function(object, condition = NULL, padj = FALSE) {
    if (is.null(assays(object)$counts)) stop("No counts found to analyze")

    if (is.null(condition)) {
        if (length(all.vars(object@design)) == 1) {
            design_var <- all.vars(object@design)[1]
            var1 <- colData(object)[[design_var]]
            cat("Welch-test on", substitute(design_var), "\n")
        }
        else {
            design_var <- colnames(colData(object))[ncol(colData(object))]
            var1 <- colData(object)[[design_var]]
            cat("Welch-test on", substitute(design_var), "\n")
        }
    }
    else {
        if (condition %in% colnames(colData(object))) {
            var1 <- colData(object)[[condition]]
            cat("Welch-test on", substitute(condition), "\n")
        }
        else (stop(paste0("Condition: ", condition, " not found in colData")))
    }
    object@results$wilcox_test_results <- run_wilcox(my_data = assays(object)$counts, var1 = var1, padj = padj)
    return(object)
})

