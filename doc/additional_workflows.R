## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(statomatic)

## ----read in norm_counts------------------------------------------------------
x <- read.csv("normalized_counts.csv")
head(x)

## ----cf norm_counts-----------------------------------------------------------
x <- cf(x)
head(x)

## ----part 1 sample_info-------------------------------------------------------
sample_info <- data.frame(samples = colnames(x), sex = c(rep("Male", 10), rep("Female", 10)), 
                          genotype = c(rep("WT", 5), rep("KO", 5), rep("WT", 5), rep("KO", 5)), 
                          group = c(rep("M_WT", 5), rep("M_KO", 5), rep("F_WT", 5), rep("F_KO", 5)))
sample_info

## ----multigroup_main 1--------------------------------------------------------
res <- multigroup_main(x, var1 = sample_info$sex, var2 = sample_info$genotype, var3 = sample_info$group)
names(res)

## ----multigroup_main 2--------------------------------------------------------
res <- multigroup_main(x, var1 = sample_info$group)
names(res)

## ----twogroup_main 1----------------------------------------------------------
res <- twogroup_main(x, var = sample_info$sex)
names(res)

## ----source code multigroup_main----------------------------------------------
multigroup_main

## ----test_norm----------------------------------------------------------------
res <- test_norm(x, sample_info$sex, sample_info$genotype)

## ----test_norm 2--------------------------------------------------------------
res <- test_norm(x, sample_info$group)

## ----run_anova----------------------------------------------------------------
anova_res <- run_anova(x, sample_info$sex, sample_info$genotype)
head(anova_res)

## ----run_tukey----------------------------------------------------------------
tukey_res <- run_tukey(x, sample_info$group)
head(tukey_res)

## ----list files---------------------------------------------------------------
#get all files in the current directory with abundance.csv in the name 
a <- list.files(".", pattern = "abundance.csv")
a

## ----read in------------------------------------------------------------------
#read them in with read.csv
read_in <- lapply(a, read.csv)
#look at first three rows of data
lapply(read_in, head, 3)

## ----cf all-------------------------------------------------------------------
#apply cf to each table
tables <- lapply(read_in, cf)

## ----make sample_info---------------------------------------------------------
#tables[[1]] is the first table in the tables list. all the tables have the same column names
sample_info <- data.frame(sample = colnames(tables[[1]]), group = c(rep("group1", 8), rep("group2", 8), rep("group3", 8)))
sample_info

## ----lapply main--------------------------------------------------------------
res <- lapply(tables, multigroup_main, sample_info$group)

## -----------------------------------------------------------------------------
#res has 6 elements, one per table we analyzed
length(res)
#each element of res has 7 elements
length(res[[1]])
names(res[[1]])

## ----results------------------------------------------------------------------
#get the first element of res (which is a list of results of the first analysis). This will be the class abundances since it is first alphabetically
x <- res[[1]]
lapply(x[1:3], head)

## ----name list----------------------------------------------------------------
#gsub .csv with nothing in list a
name_list <- gsub(".csv", "", a)
name_list

## ----collate tables-----------------------------------------------------------
#define an empty list for our new tables
collated_tables <- list()
#use a for loop to collate and save our tables
for (i in 1:length(res)) {
  #the current result were processing
  temp <- res[[i]]
  #get just the p values and fold change results 
  first_three <- temp[1:3]
  #rbind into one new table
  new_table <- rbind(first_three[[1]], first_three[[2]], first_three[[3]])
  #save the new table to the collated_tables list
  collated_tables[[i]] <- new_table
}
#set the names of collated_tables to be the names from name_list
names(collated_tables) <- name_list
#look at some of the results
lapply(collated_tables, head, 3)

## ----eval = FALSE, export-----------------------------------------------------
#  for (i in 1:length(collated_tables)) {
#    write.csv(collated_tables[[i]], paste0("results_", a[i]))
#  }

