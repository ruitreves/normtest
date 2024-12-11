## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(statomatic)

## ----data_read-in-------------------------------------------------------------
norm_counts <- read.csv("normalized_counts.csv")
head(norm_counts)

## ----column_fix---------------------------------------------------------------
norm_counts <- cf(norm_counts)
head(norm_counts)

## ----make sample info---------------------------------------------------------
sample_info <- data.frame(samples = colnames(norm_counts), sex = c(rep("Male", 10), rep("Female", 10)), 
                          genotype = c(rep("WT", 5), rep("KO", 5), rep("WT", 5), rep("KO", 5)), 
                          group = c(rep("M_WT", 5), rep("M_KO", 5), rep("F_WT", 5), rep("F_KO", 5)))
sample_info

## ----build sds part1----------------------------------------------------------
sds <- make_sds(x = norm_counts, colData = sample_info, design = ~ sex + genotype)
sds

## ----run sds_analyze----------------------------------------------------------
sds <- sds_analyze(sds)
sds

## ----all results--------------------------------------------------------------
res <- get_results(sds)
head(res)

## ----anova results------------------------------------------------------------
anova_res <- get_anova(sds)
head(anova_res)

## ----names--------------------------------------------------------------------
names(sds@results)

## ----read in------------------------------------------------------------------
x <- read.csv("two_group_normalized_counts.csv")
head(x)

## ----cf-----------------------------------------------------------------------
x <- cf(x)
head(x)

## ----sample_info--------------------------------------------------------------
sample_info <- data.frame(sample = colnames(x), group = c(rep("treated", 4), rep("untreated", 4)))
sample_info

## ----build sds part2----------------------------------------------------------
sds <- make_sds(x, colData = sample_info, design = ~ group)
sds

## ----analyze------------------------------------------------------------------
sds <- sds_analyze(sds)

## ----get results--------------------------------------------------------------
res <- get_results(sds)
head(res)

## ----@results-----------------------------------------------------------------
names(sds@results)

## -----------------------------------------------------------------------------
#read in data
x <- read.csv("species_abundance.csv")
#fix rownames
x <- cf(x)
head(x)

## -----------------------------------------------------------------------------
#make sample_info
sample_info <- data.frame(sample = colnames(x), group = c(rep("group1", 8), rep("group2", 8), rep("group3", 8)))
sample_info

## -----------------------------------------------------------------------------
#make sds
sds <- make_sds(x, colData = sample_info, design = ~ group)
sds <- sds_analyze(sds)

## -----------------------------------------------------------------------------
res <- get_results(sds)
head(res)

