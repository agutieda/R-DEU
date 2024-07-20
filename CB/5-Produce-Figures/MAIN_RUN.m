% RANDOM DISCOUNTED EXPECTED UTILITY
%
% Generate figures reported in the paper
%
% This script:
% - Generate predicted and observed choice distribution plots for the R-DEU
%   model and alternative models shown in Figure 3
% - Generate predicted and observed choice distribution plots by risk task for
%   the R-DEU model shown in Figure 4
% - Generate predicted and observed choice distribution plots by menu for the
%   R-DEU and alternative models shown in the supplementary material
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
copyfile '../0-Prepare-Data/output/' './input/'
copyfile '../1-Compute-RDEU/output/' './input/'
copyfile '..'/2-Compute-RUM/output/ './input/'
copyfile '..'/3-Compute-NLS/output/ './input/'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Figure 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure_3a_rdeu_predicted_choice_all;
figure_3b_rum_predicted_choice_all;
figure_3c_nls_predicted_choice_all;
figure_3d_rdeu_predicted_choice_menu;
figure_3e_rum_predicted_choice_menu;
figure_3f_nls_predicted_choice_menu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2) Figure 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure_4_rdeu_predicted_choice_by_task;

% Clean folders
delete('./input/*')