# R-DEU: Extensions

This folder contains the files necessary to replicate the results reported in
the section "Extensions" of the Appendix. The folders are numbered in the order
in which they should be run:

1. The folder "**0-Prepare-Data**" contains R code to load and process the file
   "dataset.dta" with the data from the replication files of the paper
   *"Estimating Time Preferences from Convex Budgets"* by James Andreoni and
   Charles Sprenger (American Economic Review, 2012).

   This Stata file can be obtained directly from the subfolder "Data" in the
   replication material of the paper available
   [here](https://www.aeaweb.org/articles?id=10.1257/aer.102.7.3333).

   With "dataset.dta" in subfolder "input", run the R script `RUN.R` to process
   the raw data and and select the sample used in the Appendix "Extensions" of
   the paper.

2. The folder "**1-Compute-Baseline**" contains Matlab code to compute the
   maximum-likelihood estimator of the joint distribution of risk and time
   preferences for the whole sample using the parametric R-DEU model described
   in Section 6.1 of the paper. Run the script `RUN.m` in the folder to compute
   these estimates.

3. The folder "**2-Compute-QMC**" contains Matlab code to compute the
   maximum-likelihood estimator of the joint distribution of risk and time
   preferences for the whole sample using the parametric R-DEU model with the
   Quasi-Montecarlo Method described in the appendix of the paper. Run the
   script `RUN.m` in the folder to compute these estimates.

4. The folder "**3-Compute-QMC-Constrained**" contains Matlab code to compute
   the constrained maximum-likelihood estimator of the joint distribution of
   risk and time preferences for the whole sample using the parametric R-DEU
   model with the Quasi-Montecarlo Method described in the appendix of the
   paper, where the distribution of the discount rate follows a truncated normal
   distribution and a Gaussian copula is used to allow for correlation between
   risk and time preferences. Run the script `RUN.m` in the folder to compute
   these estimates.

5. The folder "**4-Compute-HD**" contains Matlab code to compute the
   maximum-likelihood estimator of the joint distribution of risk and time
   preferences for the whole sample using the parametric R-DEU model extended to
   allow for present bias in discounting. Run the script `RUN.m` in the folder
   to compute these estimates.

6. The folder **5-Compute-HD-Constrained** contains Matlab code to compute the
   constrained maximum-likelihood estimator of the joint distribution of risk
   and time preferences for the whole sample using the parametric R-DEU model
   extended to allow for present bias in discounting, where the distribution of
   the discount rate follows a truncated normal distribution, the present bias
   follows a Beta distribution and a Gaussian copula is used to allow for
   correlation between preference parameters. Run the script `RUN.m` in the
   folder to compute these estimates.

7. The folder **6-Produce-Tables** reads the output from the previous steps and
   produces the Table summarizing the results and reported in the Appendix
   "Extensions".To reproduce this table, run the script `RUN.m` in the folder.

Further information can be found in the comments at the beginning of each
function and script.
