# Empirical Illustration of R-DEU Using Convex Budget Data

This folder contains the files necessary to replicate the results reported in
Tables 4 and Figures 3 and 4 in Section 6 of the paper. Its subfolders are
numbered in the order in which they should be run:

1. The subfolder **0-Prepare-Data** contains R code to load and process the file
   "dataset.dta". which contains data from the replication files of the paper
   *"Risk Preferences Are Not Time Preferences"* by James Andreoni and Charles
   Sprenger (American Economic Review, 2012). This file can be obtained directly
   from the subfolder "Data" in the replication material of the paper available,
   availale [here](https://www.aeaweb.org/articles?id=10.1257/aer.102.7.3357).

   Place the "dataset.dta" in the subsubfolder "input" and run the R script
   `RUN.R` to process the raw data and and select the sample used in Section 4
   of the paper.

2. The subfolder **1-Compute-RDEU** contains Matlab code to compute the
   maximum-likelihood estimator of of the joint distribution of risk and time
   preferences for the whole sample using the parametric R-DEU model described
   in Section 6.1 of the paper. Run the script `RUN.m` in this subfolder to
   compute these estimates and their standard errors.

3. The subfolder **2-Compute-RUM** contains Matlab code to compute the
   maximum-likelihood estimator of mean and variance of risk and time
   preferences for the whole sample using the multinomial iid-RUM model
   described in Section 6.1 of the paper. Run the script `RUN.m` in this
   subfolder to compute these estimates and their standard errors.

4. The subfolder **3-Compute-NLS** contains Matlab code to estimate risk and
   time preferences using the non-linear least squares estimator described in
   Section 6.1 of the paper. Run the script `RUN.m` in this subfolder to compute
   these estimates and their standard errors.

5. The subfolder **4-Produce-Tables** reads the output from the previous steps
   and produces Table 4 in the paper. To reproduce this table, run the script
   `RUN.m` in this subfolder after running the corresponding files in subfolders
   0 to 3.

6. The subfolder **5-Produce-Figures** reads the output from the previous steps
   and produces Figures 3 and 4 in the paper. To reproduce these figures, run
   the script `RUN.m` in this subfolder after running the corresponding files in
   subfolders 0 to 3.

Further information can be found in the comments at the beginning of each
function and script.

The code has been tested using Matlab R2024a and R version 4.4.1 on a Windows 10
machine.