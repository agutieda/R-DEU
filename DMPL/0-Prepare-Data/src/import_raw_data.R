# RANDOM DISCOUNTED EXPECTED UTILITY
# Replication Files for empirical analysis using Double-Multiple Price Lists
#
# This file:
#   - Read data from replication files of AHLR (2008)
#   - Select the variables corresponding to the lotteries and the
#     corresponding choice of each experimental subject

## 1) Load libraries and import raw dataset from the replication files ---------

rm(list = ls())
library(haven)
library(dplyr)
library(readr)
library(magrittr)

raw_data <- read_dta("./input/tmp.dta")

## 2) Select variables of interest ---------------------------------------------

selected_data <- raw_data %>%
    select(
        id,                   # ID of experimental subject
        choice,               # A dummy indicating choice in their data
        ra,                   # A dummy indicating risk task
        RowProbA,             # Probability of 1th outcome in risk tasks
        horizon,              # Delay in delay tasks
        Vprizea1,             # Payoff: 1th outcome of 1th lottery in risk tasks
        Vprizea2,             # Payoff: 2nd outcome of 1th lottery in risk tasks
        Vprizeb1,             # Payoff: 1th outcome of 2th lottery in risk tasks
        Vprizeb2,             # Payoff: 2nd outcome of 2th lottery in risk tasks
        principal,            # Current payoff in delay tasks
        DelayedP,             # Delayed payoff in delay tasks
        starts_with("Rtask"), # Risk task ID
        starts_with("DRtask") # Delay task ID
    )

## 3) Recode and summarize the variables ---------------------------------------

# We take the following convention:
#   - For risk tasks  -> L1: Safe lottery, L2: Risky Lottery
#   - For delay tasks -> L1: Early Choice, L2: Later Choice

# Given this condition, we create an indicator function "Y" that takes value 1
# when L1 is chosen 2 when L2 is chosen and 0 when indifference is reported.

# The indicator "choice" from the original dataset is coded as:
# a: early or safe choices
# b: risky or later choices
# i: indifference

selected_data %<>% mutate(
    Y = case_when(choice == "a" ~  1 ,
                  choice == "b" ~  2 ,
                  choice == "i" ~  0 ,
                  TRUE ~ NA_real_) )

# Create indicator of task set and task type
selected_data %<>% mutate(
    taskID = case_when(Rtask1  == 1 ~  1  ,
                       Rtask2  == 1 ~  2  ,
                       Rtask3  == 1 ~  3  ,
                       Rtask4  == 1 ~  4  ,
                       DRtask1 == 1 ~  5  ,
                       DRtask2 == 1 ~  6  ,
                       DRtask3 == 1 ~  7  ,
                       DRtask4 == 1 ~  8  ,
                       DRtask5 == 1 ~  9  ,
                       DRtask6 == 1 ~  10 ,
                       TRUE ~ NA_real_) )

# Create indicator of task type
selected_data %<>% mutate(
    taskType = case_when(taskID  <= 4 ~ 'Risk'   ,
                         taskID  >= 5 ~ 'Time'  ,
                         TRUE ~ NA_character_) )

# Create payoffs for each lottery consistent with our notation

scale_factor <- 1000 # Express in thousands of DKK to avoid numerical problems

selected_data %<>% mutate(
    x1_L1 = case_when(taskType == 'Risk' ~ Vprizea1  / scale_factor,
                      taskType == 'Time' ~ principal / scale_factor,
                      TRUE ~ NA_real_),
    x2_L1 = case_when(taskType == 'Risk' ~ Vprizea2  / scale_factor,
                      taskType == 'Time' ~ principal / scale_factor,
                      TRUE ~ NA_real_),
    x1_L2 = case_when(taskType == 'Risk' ~ Vprizeb1  / scale_factor,
                      taskType == 'Time' ~ DelayedP  / scale_factor,
                      TRUE ~ NA_real_),
    x2_L2 = case_when(taskType == 'Risk' ~ Vprizeb2  / scale_factor,
                      taskType == 'Time' ~ DelayedP  / scale_factor,
                      TRUE ~ NA_real_) )

# Create probabilities and delay times consistent with our notation
selected_data %<>% mutate(
    p_L1 = case_when(taskType == 'Risk' ~ RowProbA ,
                     taskType == 'Time' ~ 1 ,
                     TRUE ~ NA_real_),
    p_L2 = case_when(taskType == 'Risk' ~ RowProbA ,
                     taskType == 'Time' ~ 1 ,
                     TRUE ~ NA_real_),
    t_L1 = case_when(taskType == 'Risk' ~ 0  ,
                     taskType == 'Time' ~ 30/360 ,
                     TRUE ~ NA_real_),
    t_L2 = case_when(taskType == 'Risk' ~ 0  ,
                     taskType == 'Time' ~ (30+horizon)/360 ,
                     TRUE ~ NA_real_) )

# Notice that in their data, p1_L1 = p1_L2 and that the payoff in the
# early choice is always delivered 30 days in the future

## 4) Keep variables of interest -----------------------------------------------

Data_AHLR <- selected_data %>%
    select(id, x1_L1, x2_L1, x1_L2, x2_L2, p_L1, p_L2, t_L1, t_L2, Y,
           taskID, taskType)

## 5) Create menu table --------------------------------------------------------

# Add unique menu ID
Data_AHLR %<>%
    group_by(x1_L1, x2_L1, x1_L2, x2_L2, p_L1, p_L2, t_L1, t_L2) %>%
    mutate(menuID = cur_group_id()) %>%
    ungroup()

# Identify unique menus (make sure there are 100 unique menus)
menuTab <- Data_AHLR %>%
    select(menuID, taskID, taskType,
           x1_L1, x2_L1, x1_L2, x2_L2, p_L1, p_L2, t_L1, t_L2 ) %>%
    distinct()

# Sort by taskID, probability and payoff and create new ID
menuTab %<>%
    arrange(taskID,p_L1,x2_L2,t_L2) %>%
    mutate(new_menuID = row_number())

# Store this new ID to match to observations
newID <- menuTab %>% select(menuID,new_menuID)

# Drop from menuTab
menuTab %<>% select(-menuID) %>%  rename(menuID = new_menuID)


## 6) Create obs and menu tables -----------------------------------------------

# Attach to original tab and discard old menuID
Data_AHLR <- left_join(Data_AHLR, newID, by = "menuID") %>%
    select(-menuID) %>%
    rename(menuID = new_menuID)

# Extract information about observations
obsTab <- Data_AHLR %>%
    select(subjectID = id, menuID, Y ) %>%
    arrange(subjectID, menuID)

## 7) Export -------------------------------------------------------------------
write_csv(menuTab, "./temp/rawMenuTab.csv")
write_csv(obsTab , "./temp/rawObsTab.csv")
