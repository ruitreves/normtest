# normtest
This package was developed to address the need for a more robust statistical framework to process RNA sequencing data. 

## Introduction
Certain experiemental designs are often suited to the preplanned statistical analysis that will take place after data is collected. It is vital, however, to apply the appropriate statistical techniques in order to ensure the validity of the results. The purpose of this package is to provide a high throughput way to sort data based on which tests are appropriate, and carry out those tests. 

As an example, consider a researcher who hopes to perform a two way anova test on a dataset. The researcher may plug their data in and generate statistical results, but without consideration for the mathematical assumptions of the two way anova test, their results could be erroneous. To be suitable for a two way anova, the residuals of the data must be normally distributed and the data must be homoscedastic. The normtest package is equipped to address these assumptions and sort the data based on which tests can be appropriately applied. 

## Details
The main functions, "rui" and "two_group_rui", pull together the functions of this package into the appropriate order, but the functions are available to use individually. The general process is as follows:

1. For each row in a data object, a linear model is created based on the specified factors (groups, etc.). From this linear model, the residuals are tested using a Shapiro Normality Test with alpha = 0.05. Rows that are deemed by this criteria to have normally distributed residuals are then tested by Levene's Test to assess the equality of variances (homoscedacity) amongst the specified factors with alpha = 0.05. Rows that are deemed by this criteria to be homoscedastic are used to construct a new data object, "tier1". Rows that are deemed to be heteroscedastic (non-homoscedastic) are used to construct the data object "tier2". Rows that are deemed to be non-normally distributed are used to construct the data object "tier3"
2. In the case of multiple groups (> 2 groups):
     Tier1 is analyzed by a two way anova test, reporting a p value for factor1, factor2, and the interaction of factor1 and factor2. A Tukey multi comparison test is done based on factor3, which is assumed to be some linear combination of factors 1 and 2 (e.g., if factor1 = A, factor2 = B, then factor 3 should be something like AB, either literally or symbolically). P values are reported for each groupwise comparison.
     Tier2 is analyzed by a Welch's t-test based on factor3, and a p value is reported for this test. A Dunnett T3 multicomparison test is done based also on factor3, and a p value is reported for each groupwise comparison. 
     Tier3 is analyzed by a Kruskal-Wallis rank sums test and a p value is reported. This is followed by a Dunn multicomparison test and a p value is reported for each groupwise comparison.
   In the case of two groups:
     Tier1 is analyzed by a Student's t-test and a p value is reported.
     Tier2 is analyzed by a Welch's t-test and a p value is reported.
     Tier3 is analyzed by a Wilcoxon/Mann-Whitney test and a p value is reported. 
3. After the statisitcal tests are performed, log2 fold changes are calculated. Note that the program will sort the group names alphabetically before calculating the fold changes. The fold changes will be reported as log2FoldChange_GroupA_vs_GroupB, which indicates that the following operation was performed: 
log2(mean(groupA) / mean(groupB))
4. The program will then order the statstical groupwise comparison columns to be next to their fold change columns, and write csv files which will specify which tests were performed. Also note that global variables "tier1", "tier2", and "tier3" will be assigned to the global R environment. This is so you can view them if you choose, but be aware that this will overwrite any other variables with those names in your environment. 
    