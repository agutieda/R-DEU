# RANDOM DISCOUNTED EXPECTED UTILITY
# Replication Files for empirical analysis using convex budgets
#
# This file:
#   - Run all scripts that prepare the dataset used in empirical analysis of
#     data from convex time budgets (CTB) designs. We use data from the
#     experiment in:
#    "Risk Preferences Are Not Time Preferences"
#     by James Andreoni and Charles Sprenger. 
#     American Economic Review, 2012.
#     https://www.aeaweb.org/articles?id=10.1257/aer.102.7.3357
#
# Input:
#   - dataset.dta: Stata file with data used in the original paper. This file
#     can be obtained directly from the subfolder "Data" in replication material
#     available at the link above.
#
# Output:
#   - menuTab.csv: Table describing each menu available in the dataset
#   - obsTab.csv: Dataset of observed choice by each subject in each menu of
#     the risk and time tasks in the dataset
#
# Written by Jose Apesteguia, Miguel A. Ballester, and Angelo Gutierrez-Daza
# July 2024
# Tested using R.4.4.1

rm(list = ls())

## 0) Create "temp" and "output" folders for intermediate inputs and output ----
unlink("./output", recursive = TRUE)
dir.create("./output")

## 1) Load libraries and import raw dataset from the replication files ---------
rm(list = ls())
library(haven)
library(dplyr)
library(readr)
library(magrittr)

raw_data <- read_dta("input/dataset.dta")

## 2) Select variables of interest ---------------------------------------------

selected_data <- raw_data %>%
    select(
        subjectID = labnumber,   # ID of experimental subject
        menuID    = budgetnum,   # ID of menu faced
        blockID   = choiceblock, # ID of block faced
        tok_soon  = v,           # Tokens allocated to sooner payment
        q_t       = rate1,       # Sooner payment token rate
        p_t       = perc1,       # Probability (in %) of receiving sooner payment
        p_tk      = perc2,       # Probability (in %) of receiving later payment
        k_days    = k,           # Delay length in days
    )

## 3) Add additional variables describing each menu ----------------------------

# Add indicator of risk condition
selected_data %<>% mutate(
    taskID = case_when((blockID  == 1) | (blockID  == 2)  ~  1  ,
                       (blockID  == 3) | (blockID  == 4)  ~  2  ,
                       (blockID  == 5) | (blockID  == 6)  ~  3  ,
                       (blockID  == 7) | (blockID  == 8)  ~  4  ,
                       (blockID  == 9) | (blockID  == 10)  ~  5  ,
                       (blockID  == 11) | (blockID  == 12)  ~  6  ,
                       TRUE ~ NA_real_
    ))

# Add token rate in t+k (in this experiment, it's always 0.2)
selected_data %<>% mutate(q_tk = 0.2)

# Add token budget (in this experiment, it's always 100)
selected_data %<>% mutate(tok_budget = 100)


# Write probabilities in decimal form and the probability ratio
selected_data %<>%  mutate(p_t  = p_t / 100 ,
                           p_tk = p_tk / 100 ,
                           p_ratio = p_tk / p_t)

# Add payoff dates in years
selected_data %<>% mutate(k    = k_days / 360 ,
                          t_t  = 7 / 360 ,
                          t_tk = (7 + k_days) / 360)

# Add implied return rate and implied annual rate
selected_data %<>% mutate(rr   = log(q_tk / q_t),
                          rr_y = ( (q_tk / q_t) ^ (360 / k_days) - 1)*100 )

# Compute implicit consumption choice in t and t+k
selected_data %<>% mutate(c_t  = q_t  * tok_soon,
                          c_tk = q_tk * (tok_budget - tok_soon))

# Compute share of tokens allocated to the future
# (alpha = 1 -> all tokens in t+k)
selected_data %<>% mutate(alpha_original = (1 - tok_soon / tok_budget))

# Compute alpha approximated to closest multiple of 0.05
selected_data %<>% mutate(
    alpha = case_when(
        (alpha_original >= 0)     & (alpha_original <= 0.025) ~ 0.00,
        (alpha_original  > 0.025) & (alpha_original <= 0.075) ~ 0.05,
        (alpha_original  > 0.075) & (alpha_original <= 0.125) ~ 0.10,
        (alpha_original  > 0.125) & (alpha_original <= 0.175) ~ 0.15,
        (alpha_original  > 0.175) & (alpha_original <= 0.225) ~ 0.20,
        (alpha_original  > 0.225) & (alpha_original <= 0.275) ~ 0.25,
        (alpha_original  > 0.275) & (alpha_original <= 0.325) ~ 0.30,
        (alpha_original  > 0.325) & (alpha_original <= 0.375) ~ 0.35,
        (alpha_original  > 0.375) & (alpha_original <= 0.425) ~ 0.40,
        (alpha_original  > 0.425) & (alpha_original <= 0.475) ~ 0.45,
        (alpha_original  > 0.475) & (alpha_original <= 0.525) ~ 0.50,
        (alpha_original  > 0.525) & (alpha_original <= 0.575) ~ 0.55,
        (alpha_original  > 0.575) & (alpha_original <= 0.625) ~ 0.60,
        (alpha_original  > 0.625) & (alpha_original <= 0.675) ~ 0.65,
        (alpha_original  > 0.675) & (alpha_original <= 0.725) ~ 0.70,
        (alpha_original  > 0.725) & (alpha_original <= 0.775) ~ 0.75,
        (alpha_original  > 0.775) & (alpha_original <= 0.825) ~ 0.80,
        (alpha_original  > 0.825) & (alpha_original <= 0.875) ~ 0.85,
        (alpha_original  > 0.875) & (alpha_original <= 0.925) ~ 0.90,
        (alpha_original  > 0.925) & (alpha_original <= 0.975) ~ 0.95,
        (alpha_original  > 0.975) & (alpha_original <= 1.000) ~ 1.00,
        TRUE ~ NA_real_
    )
)

# Add indicator function between 1 and 21 for the corresponding alpha
selected_data %<>% mutate(
    Y  =  case_when(
        alpha == 0.00 ~ 1,
        alpha == 0.05 ~ 2,
        alpha == 0.10 ~ 3,
        alpha == 0.15 ~ 4,
        alpha == 0.20 ~ 5,
        alpha == 0.25 ~ 6,
        alpha == 0.30 ~ 7,
        alpha == 0.35 ~ 8,
        alpha == 0.40 ~ 9,
        alpha == 0.45 ~ 10,
        alpha == 0.50 ~ 11,
        alpha == 0.55 ~ 12,
        alpha == 0.60 ~ 13,
        alpha == 0.65 ~ 14,
        alpha == 0.70 ~ 15,
        alpha == 0.75 ~ 16,
        alpha == 0.80 ~ 17,
        alpha == 0.85 ~ 18,
        alpha == 0.90 ~ 19,
        alpha == 0.95 ~ 20,
        alpha == 1.00 ~ 21,
        TRUE ~ NA_real_
    )
)

## 4) Keep variables of interest -----------------------------------------------
Data_AS <- selected_data %>%
    select(
        subjectID,
        menuID,
        blockID,
        taskID,
        Y,
        alpha,
        alpha_original,
        c_t,
        c_tk,
        t_t,
        t_tk,
        p_t,
        p_tk,
        q_t,
        q_tk,
        k,
        p_ratio,
        rr,
        rr_y
    ) %>%
    arrange(subjectID, menuID)


## 5) Create menu table --------------------------------------------------------

# Identify unique menus (make sure there are 84 unique menus)
menuTab_AS <- Data_AS %>%
    select(menuID,
           blockID,
           taskID,
           t_t,
           t_tk,
           p_t,
           p_tk,
           q_t,
           q_tk,
           k,
           p_ratio,
           rr,
           rr_y) %>%
    distinct() %>%
    arrange(menuID)

## 6) Create obs table ---------------------------------------------------------
obsTab_AS <- Data_AS %>%
    select(subjectID,
           menuID,
           blockID,
           taskID,
           Y,
           alpha,
           alpha_original,
           c_t,
           c_tk) %>%
    arrange(subjectID, menuID)

## 7) Create table with subject info -------------------------------------------
subjectTab_AS <- obsTab_AS %>%
    select(subjectID,alpha) %>%
    group_by(subjectID) %>%
    summarize(
        nAllPresent = sum(alpha==0),
        nAllFuture  = sum(alpha==1),
        nCorner     = nAllPresent + nAllFuture,
        nInterior   = sum(alpha>0 & alpha<1),
        nTotal      = nCorner + nInterior
    )

## 7) Export -------------------------------------------------------------------
write_csv(menuTab_AS , "./output/menuTab.csv")
write_csv(obsTab_AS  , "./output/obsTab.csv")
write_csv(subjectTab_AS  , "./output/subjectTab.csv")
