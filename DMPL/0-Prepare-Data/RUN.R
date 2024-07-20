# RANDOM DISCOUNTED EXPECTED UTILITY
# Replication Files for empirical analysis using Double-Multiple Price Lists
#
# This file:
#   - Run all scripts that prepare the dataset used in empirical analysis of
#     double-multiple price list (DMPL) designs. We use a subsample from the
#     experiment in:
#     "Eliciting Risk And Time Preferences"
#     by Steffen Andersen, Glenn Harrison, Morten Lau, and Elisabet Rutstroml
#     Econometrica, 2008. Available at:
#     https://www.econometricsociety.org/publications/econometrica/2008/05/01/eliciting-risk-and-time-preferences
#
# Input:
#   - tmp.dta: Stata file with data used in the original paper.
#     This can be obtained by running script "RiskAndTime.do" from the 
#     replication material available in the link above, while commenting lines
#     37-42 of that script to prevent file deletion
#
# Output:
#   - menuTab.csv: Table describing each menu available in the dataset
#   - obsTab.csv: Dataset of observed choice by each subject in each menu of
#     the risk and time tasks in the dataset
#   - Select_Sample.html: HTML file describing in detail each step of the
#     sample selection process
#
# Written by Jose Apesteguia, Miguel A. Ballester, and Angelo Gutierrez-Daza
# July 2024
# Tested using R.4.4.1

rm(list = ls())

## 0) Create "temp" and "output" folders for intermediate inputs and output ----
unlink("./temp", recursive = TRUE)
dir.create("./temp")
unlink("./output", recursive = TRUE)
dir.create("./output")

## 1) Load libraries and import raw dataset from the replication files ---------
source("./src/import_raw_data.R")

## 2) Select sample used for estimation ----------------------------------------
quarto::quarto_render("./src/select_sample.qmd")

## 3) Clean up -----------------------------------------------------------------
file.rename("./src/select_sample.html", "./output/select_sample.html")
unlink("./temp", recursive = TRUE)
