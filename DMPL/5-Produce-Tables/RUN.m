% RANDOM DISCOUNTED EXPECTED UTILITY
%
% Generate tables reported in the paper
%
% This script:
% - Generate Table 1 reporting population-level estimates of risk and time 
%   preferences uisng the R-DEU and alternative models with DMPL data.
% - Generate Table 2 reporting summary statistics of individual-level estimates 
%   of mean risk aversion and discounting using the R-DEU and alternative models
%   with DMPL data.
% - Generate Table 3 reporting summary statistics of individual-level estimates 
%   of std. dev. of risk aversion and discounting using the R-DEU and alternative
%   models with DMPL data.
% - Generate first two columns of Appendix Table reporting population-level
%   estimates of risk and time preferences using the R-DEU model under different
%   levels of background consumption.
%
% Written by Jose Apesteguia, Miguel A. Ballester, and Angelo Gutierrez-Daza
% July 2024
% Tested using Matlab 2024a

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0) Initialize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;
addpath([pwd,'/src']);

% Clean folders
delete('./input/*')
delete('./output/*')

% Import data from corresponding folder
copyfile '../1-Compute-SPE/output/' './input/'
copyfile '../2-Compute-RDEU/output/' './input/'
copyfile '../3-Compute-LUCE/output/' './input/'
copyfile '../4-Compute-WILCOX/output/' './input/'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Table 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table_1_mle_population;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2) Table 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table_2_mle_subject_mu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3) Table 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table_3_mle_subject_sigma;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4) Column 1 and 2 of Appendix Table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table_Appendix_mle_omega_c12;

% Clean folders
delete('./input/*')
