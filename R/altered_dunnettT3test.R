#'@export 
AdunnettT3Test<- function(x, g, ...) {
    #this function was altered from its original form on 12/2/24 by Rui Treves
    #The changes are as follows:
    # 1. on line 47, within the compare.stats function, "if (is.nan(tval)) {tval <- 0}" was added.
    # 2. on line 58, within the getDF function, "if (is.nan(df)) {df <- 0}" was added
    #These changed allow for two groups to have all the same values and the function will not 
    #throw an error, but instead return a p-value of 1 for such a group comparison. 
        ## taken from stats::kruskal.test

    if (is.list(x)) {
        if (length(x) < 2L)
            stop("'x' must be a list with at least 2 elements")
        DNAME <- deparse(substitute(x))
        x <- lapply(x, function(u) u <- u[complete.cases(u)])
        k <- length(x)
        l <- sapply(x, "length")
        if (any(l == 0))
            stop("all groups must contain data")
        g <- factor(rep(1 : k, l))
        x <- unlist(x)
    }
    else {
        if (length(x) != length(g))
            stop("'x' and 'g' must have the same length")
        DNAME <- paste(deparse(substitute(x)), "and",
                       deparse(substitute(g)))
        OK <- complete.cases(x, g)
        x <- x[OK]
        g <- g[OK]
        if (!all(is.finite(g)))
            stop("all group levels must be finite")
        g <- factor(g)
        k <- nlevels(g)
        if (k < 2)
            stop("all observations are in the same group")
    }

    ## prepare dunnettT3 test
    ni <- tapply(x, g, length)
    n <- sum(ni)
    xi <- tapply(x, g, mean)
    s2i <- tapply(x, g, var)

    s2in <- 1 / (n - k) * sum(s2i * (ni - 1))

    compare.stats <- function(i,j) {
        dif <- xi[i] - xi[j]
        A <- (s2i[i] / ni[i] + s2i[j] / ni[j])
        tval <- dif / sqrt(A)
        if (is.nan(tval)) {tval <- 0}
        return(tval)
    }

    PSTAT <- pairwise.table(compare.stats,levels(g), p.adjust.method="none" )

    getDf <- function(i, j){
        A <- (s2i[i] / ni[i] + s2i[j] / ni[j])
        df <- A^2 / (s2i[i]^2 / (ni[i]^2 * (ni[i] - 1)) +
                    s2i[j]^2 / (ni[j]^2 * (ni[j] - 1)))
        if (is.nan(df)) {df <- 0}
        return(df)
        }

    DF <- pairwise.table(getDf, levels(g), p.adjust.method="none" )

    ## prepare matrix
    m <- k * (k - 1) / 2 # number of comparisons
    cr <- diag(m)
    pval <- rep(NA, m)

    ## use Studentized Maximum Modulus Distribution
    ## equals Two-sided Multivariate t distribution
    ## use mvtnorm critical values

    df <- as.vector(DF)
    df <- df[!is.na(df)]
    df <- round(df, digits = 0) # round to integer
    pstat <- as.vector(PSTAT)
    pstat[pstat == NaN] <- 0
    pstat <- pstat[!is.na(pstat)]

    for (i in 1:m){
        lo <- -rep(abs(pstat[i]), m)
        up <- rep(abs(pstat[i]), m)
        pval[i] <- 1 - mvtnorm::pmvt(lower = lo,
                            upper = up,
                            df = df[i],
                            corr = cr)
    }

    ## create output matrix
    PVAL <- PSTAT
    PVAL[!is.na(PVAL)] <- pval

    DIST <- "t"
    alternative <- "two.sided"
    METHOD <- paste0("Dunnett's T3 test for multiple comparisons\n",
                     "\t\twith unequal variances")
    p.adjust.method <- "single-step"
    MODEL <- data.frame(x, g)
    ans <- list(method = METHOD, data.name = DNAME, p.value = PVAL,
                statistic = PSTAT, p.adjust.method = p.adjust.method,
                model = MODEL, dist = DIST, alternative = alternative)
    class(ans) <- "PMCMR"
    ans
}
