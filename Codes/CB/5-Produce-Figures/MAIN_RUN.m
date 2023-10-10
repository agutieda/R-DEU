% RANDOM DISCOUNTED EXPECTED UTILITY
%
% Generate figures reported in the paper
%
% This script:
% - Generate predicted and observed choice distribution plots for the RDEU model
%   and alternative models shown in Figure 2
% - Generate predicted and observed choice distribution plots by risk task for
%   the RDEU model shown in Figure 4
% - Generate predicted and observed choice distribution plots by menu for the
%   RDEU and alternative models shown in the supplementary material
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

% Clean folders
delete('./input/*')
delete('./output/*')

% Import data from corresponding folder
copyfile '../0-Prepare-Data/output/' './input/'
copyfile '../1-Compute-RDEU/output/' './input/'
copyfile '..'/2-Compute-RUM/output/ './input/'
copyfile '..'/3-Compute-NLS/output/ './input/'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Figure 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure_2a_rdeu_predicted_choice_all;
figure_2b_rum_predicted_choice_all;
figure_2c_nls_predicted_choice_all;
figure_2d_rdeu_predicted_choice_menu;
figure_2e_rum_predicted_choice_menu;
figure_2f_nls_predicted_choice_menu;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2) Figure 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure_4_rdeu_predicted_choice_by_task;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3) Supplementary Material
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure_5a_rdeu_predicted_choice_all_menus;
figure_5b_rum_predicted_choice_all_menus;
figure_5c_nls_predicted_choice_all_menus;

% Clean folders
delete('./input/*')