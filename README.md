# statomatic
This package was developed to provide a straight-forward and high-throughput method for analyzing RNA-sequencing data, but is applicable to other numeric data.

## Introduction
Certain experimental designs are often suited to the preplanned statistical analysis that will take place after data is collected. It is vital, however, to apply the appropriate statistical techniques in order to ensure the validity of the results. The purpose of this package is to provide a high-throughput way to sort data based on which tests are appropriate, and carry out those tests. 

As an example, consider a researcher who hopes to perform a two way anova test on a dataset. The researcher may plug their data in and generate statistical results, but without consideration for the mathematical assumptions of the two way anova test, their results could be erroneous. To be suitable for a two way anova, the residuals of the data must be normally distributed and the data must be homoscedastic. The statomatic package is equipped to address these assumptions and sort the data based on which tests can be appropriately applied. 

## Details
The main functions, "multigroup_main" and "twogroup_main", pull together the functions of this package into the appropriate order, but the functions are available to use individually. The general process is as follows:

1. For each row in a data object, a generalized linear model is created based on the specified factors (groups, etc.). From this model, the residuals are tested using a Shapiro Normality Test with alpha = 0.05. Rows that are deemed by this criteria to have normally distributed residuals are then tested by Levene's Test to assess the equality of variances (homoscedacity) amongst the specified factors with alpha = 0.05. Rows that are deemed by this criteria to be homoscedastic are used to construct a new data object, "ne" (for normal, equal variance). Rows that are deemed to be heteroscedastic (non-homoscedastic) are used to construct the data object "nu" (for normal, unequal variance). Rows that are deemed to be non-normally distributed are used to construct the data object "nn" (for non-normal).
2. In the case of multiple groups (> 2 groups):
     ne is analyzed by a two way anova test, reporting a p value for factor1, factor2, and the interaction of factor1 and factor2 if applicable. A Tukey multi comparison test is done based on factor3, which is assumed to be some linear combination of factors 1 and 2 (e.g., if factor1 = A, factor2 = B, then factor 3 should be something like AB, either literally or symbolically). P values are reported for each groupwise comparison.
     nu is analyzed by a Welch's t-test based on factor3, and a p value is reported for this test. A Dunnett T3 multicomparison test is done based also on factor3, and a p value is reported for each groupwise comparison. 
     nn is analyzed by a Kruskal-Wallis rank sums test and a p value is reported. This is followed by a Dunn multicomparison test and a p value is reported for each groupwise comparison.
   In the case of two groups:
     ne is analyzed by a Student's t-test and a p value is reported.
     nu is analyzed by a Welch's t-test and a p value is reported.
     nn is analyzed by a Wilcoxon/Mann-Whitney test and a p value is reported. 
3. After the statistical tests are performed, log2 fold changes are calculated. Note that the program will sort the group names alphabetically before calculating the fold changes. The fold changes will be reported as log2FoldChange_GroupA_vs_GroupB, which indicates that the following operation was performed: log2(mean(groupA) / mean(groupB))
4. Statomatic will then order the statistical groupwise comparison columns to be next to their fold change columns.


    