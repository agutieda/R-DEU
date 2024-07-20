# R-DEU: Double Multiple Price Lists

This folder contains the files necessary to replicate the results reported in
Tables 1-3 and Figures 1 and 2 in the paper. The folders are numbered in the
order in which they should be run:

1. The folder "**0-Prepare-Data**" contains R code to load and process the file
   "tmp.dta" containing data from the replication files of the paper *"Eliciting
   Risk And Time Preferences"* by Steffen Andersen, Glenn Harrison, Morten Lau,
   and Elisabet Rutstroml (Econometrica, 2008).

   To obtain this file, download the replication material of the paper,
   available
   [here](https://www.econometricsociety.org/publications/econometrica/2008/05/01/eliciting-risk-and-time-preferences),
   and run the script "RiskAndTime.do" while commenting lines 37-42 of the
   script to prevent file deletion. The script will create the desired file
   "tmp.dta".

   With "tmp.dta" in subfolder "input", run the R script `RUN.R` to process the
   raw data and and select the sample used in Section 4 of the paper.

2. The folder "**1-Compute-SPE**" contains Matlab code to compute the
   semi-parametric estimator (SPE) of the mean and variance of risk and time
   preferences described in Section 4.1 of the paper. Run the script `RUN.m` in
   the folder to compute these estimates.

3. The folder "**2-Compute-RDEU**" contains Matlab code to compute the
   maximum-likelihood estimator of of the joint distribution of risk and time
   preferences using the parametric R-DEU model described in Section 4.2 of the
   paper. The estimation is repeated for the whole sample and for each
   individual using multiple initial values. Run the script `RUN.m` in the
   folder to compute these estimates.

4. The folder "**3-Compute-LUCE**" contains Matlab code to compute the
   maximum-likelihood estimator of of the mean and variance of risk and time
   preferences using the LUCE model described in Section 4.2 of the paper. The
   estimation is repeated for the whole sample and for each individual using
   multiple initial values. Run the script `RUN.m` in the folder to compute
   these estimates.

5. The folder "**4-Compute-WILCOX**" contains Matlab code to compute the
   maximum-likelihood estimator of of the mean and variance of risk and time
   preferences using the WILCOX model described in Section 4.2 of the paper.The
   estimation is repeated for the whole sample and for each individual using
   multiple initial values. Run the script `RUN.m` in the folder to compute
   these estimates.

6. The folder "**5-Produce-Tables**" reads the output from the previous steps
   and produces Tables 1, 2, and 3 in the paper. To reproduce them, run the
   script `RUN.m` in the folder.

7. The folder "**6-Produce-Figures**" reads the output from the previous steps
   and produces Figures 1 and 2 in the paper. To reproduce them, run the script
   `RUN.m` in the folder.

Further information can be found in the comments at the beginning of each
function and script.
