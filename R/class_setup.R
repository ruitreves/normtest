#' @import SummarizedExperiment
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' 
setClass(
    "StatomaticDataset",
    contains = "SummarizedExperiment",
    slots = list(
        results = "list",
        design = "formula"
    )
)

#' Create a StatomaticDataset (sds)
#' @param x a numeric matrix or data.frame
#' @param colData a data.frame of mappings for columns of x to experimental conditions
#' @param rowData row metadata
#' @param design a design formula containing column names from colData
#' @return a StatomaticDataset
#' @export
make_sds <- function(x, colData = NULL, rowData = NULL, design = ~ 1) {
    if (!is.matrix(x)) {
        x <- as.matrix(x)
        if (!is.matrix(x)) { stop("Counts must be a matix") }
    }
    if (is.null(colData)) stop("colData is required to build the StatomaticDataset")
    if (design == ~1) warning("Design is ~1... it is recommended to define at least a one factor design")
    
    design_vars = all.vars(design)
    if (!all(design_vars %in% colnames(colData))) {
        stop("The design formula contains variables not found in the column data")
    }

    se <- SummarizedExperiment(assays = list(counts = x), rowData = rowData, colData = colData)
    new("StatomaticDataset", se, design = design, results = list())
}

#add accessors 
setGeneric("get_counts", function(object) standardGeneric("get_counts"))
#' @export
setMethod("get_counts", "StatomaticDataset", function(object) {
    if (is.null(assays(sds)$counts)) { stop("No counts found") }
    return(data.frame(assays(sds)$counts))
})

setGeneric("get_ne", function(object) standardGeneric("get_ne"))
#' @export
setMethod("get_ne", "StatomaticDataset", function(object) {
    if (is.null(object@results$ne)) { stop("Result not found, have you analyzed data yet?")}
    return(object@results$ne)
})
setGeneric("get_nu", function(object) standardGeneric("get_nu"))
#' @export
setMethod("get_nu", "StatomaticDataset", function(object) {
    if (is.null(object@results$nu)) { stop("Result not found, have you analyzed data yet?")}
    return(object@results$nu)
})
setGeneric("get_nn", function(object) standardGeneric("get_nn"))
#' @export
setMethod("get_nn", "StatomaticDataset", function(object) {
    if (is.null(object@results$nn)) { stop("Result not found, have you analyzed data yet?")}
    return(object@results$nn)
})
setGeneric("get_t_test", function(object) standardGeneric("get_t_test"))
#' @export
setMethod("get_t_test", "StatomaticDataset", function(object) {
    if(is.null(object@results$t_test_results)) {
        stop("Data has not yet been analyzed")
    }
    return(object@results$t_test_results)
})
setGeneric("get_welch_test", function(object) standardGeneric("get_welch_test"))
#' @export
setMethod("get_welch_test", "StatomaticDataset", function(object) {
    if(is.null(object@results$welch_test_results)) {
        stop("Data has not yet been analyzed")
    }
    return(object@results$welch_test_results)
})
setGeneric("get_wilcox_test", function(object) standardGeneric("get_wilcox_test"))
#' @export
setMethod("get_wilcox_test", "StatomaticDataset", function(object) {
    if(is.null(object@results$wilcox_test_results)) {
        stop("Data has not yet been analyzed")
    }
    return(object@results$wilcox_test_results)
})
setGeneric("get_results", function(object) standardGeneric("get_results"))
#' @export
setMethod("get_results", "StatomaticDataset", function(object) {
    if(is.null(object@results$all_results)) {
        stop("Data has not yet been analyzed")
    }
    return(object@results$all_results)
})
setGeneric("get_anova", function(object) standardGeneric("get_anova"))
#' @export
setMethod("get_anova", "StatomaticDataset", function(object) {
    if(is.null(object@results$anova_results)) {
        stop("Data has not yet been analyzed")
    }
    return(object@results$anova_results)
})
setGeneric("get_kw", function(object) standardGeneric("get_kw"))
#' @export
setMethod("get_kw", "StatomaticDataset", function(object) {
    if(is.null(object@results$kw_results)) {
        stop("Data has not yet been analyzed")
    }
    return(object@results$kw_results)
})
setGeneric("get_foldchange", function(object) standardGeneric("get_foldchange"))
#' @export
setMethod("get_foldchange", "StatomaticDataset", function(object) {
    if(is.null(object@results$foldchange)) {
        stop("No fold changes found")
    }
    return(object@results$foldchange)
})
setGeneric("get_tukey_test", function(object) standardGeneric("get_tukey_test"))
#' @export
setMethod("get_tukey_test", "StatomaticDataset", function(object) {
    if(is.null(object@results$tukey_test_results)) stop("Tukey results not found")
    return(object@results$tukey_test_results)
})
setGeneric("get_dunnett_test", function(object) standardGeneric("get_dunnett_test"))
#' @export
setMethod("get_dunnett_test", "StatomaticDataset", function(object) {
    if(is.null(object@results$dunnett_test_results)) stop("Dunnett results not found")
    return(object@results$dunnett_test_results)
})
setGeneric("get_dunn_test", function(object) standardGeneric("get_dunn_test"))
#' @export
setMethod("get_dunn_test", "StatomaticDataset", function(object) {
    if(is.null(object@results$dunn_test_results)) stop("Dunn results not found")
    return(object@results$dunn_test_results)
})