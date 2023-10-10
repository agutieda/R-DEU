% RANDOM DISCOUNTED EXPECTED UTILITY
%
% Generate tables reported in the paper
%
% This script:
% - Generate table 3 reporting ML estimates of the RDEU and alternative models
%   using CB data and pooling all observations in the sample, assuming a 
%   baseline level of background consumption of 5 USD
% - Generate columns 3 and 4 of Table 5 reporting ML estimates of the RDEU model
%   using CB data for two different levels of background consumption
%
% Written by Jose Apesteguia, Miguel A. Ballester, and Angelo Gutierrez-Daza
% October 2023
% Tested using Matlab 2023a

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
copyfile '../1-Compute-RDEU/output/' './input/'
copyfile '..'/2-Compute-RUM/output/ './input/'
copyfile '..'/3-Compute-NLS/output/ './input/'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Table 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table_3_mle_cb;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2) Columns 3 and 4 of Table 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

table_5_mle_omega_c34;

% Clean folders
delete('./input/*')