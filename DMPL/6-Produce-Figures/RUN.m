% RANDOM DISCOUNTED EXPECTED UTILITY
%
% Generate DMPL figures reported in the paper
%
% This script:
% - Generate scatterplots of estimates of mu_r and mu_delta by subject shown in
%   Figure 1
% - Generate scatterplots of estimates of sigma_r and sigma_delta by subject
%   shown in Figure 2
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

% Clean folders
delete('./input/*')
delete('./output/*')

% Import data from corresponding folder
copyfile '../1-Compute-SPE/output/' './input/'
copyfile '../2-Compute-RDEU/output/' './input/'
copyfile '../3-Compute-LUCE/output/' './input/'
copyfile '../4-Compute-WILCOX/output/' './input/'

rdeuTab = readtable('./input/rdeu_subject_mle.csv');
rumLTab = readtable('./input/luce_subject_mle.csv');
rumWTab = readtable('./input/wilcox_subject_mle.csv');
spRTab = readtable('./input/spe_subject_ra.csv','ReadRowNames',true);
spTTab = readtable('./input/spe_subject_dr.csv','ReadRowNames',true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Fig 1: Scatter of mu_r and mu_delta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure_1_scatter_mu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Fig 2: Scatter of sigma_r and sigma_delta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure_2_scatter_sigma;


% Clean folders
delete('./input/*')