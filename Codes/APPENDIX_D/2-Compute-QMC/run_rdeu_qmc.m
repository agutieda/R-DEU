% RANDOM DISCOUNTED EXPECTED UTILITY
%
% Estimate RDEU model using Convex Budget (CB) data
%
% This script:
% - Estimate risk and time preferences at population level
% - Use CB data from:
%  "Estimating Time Preferences from Convex Budgets"
%   by James Andreoni and Charles Sprenger.
%   American Economic Review, (2012)
% - Assume RDEU model where (ra,dr) follow a bivariate normal distribution
% - Use quasi-monte carlo methods
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
rng(1);
addpath([pwd,'/src']);

% Clean folders
delete('./input/*')
delete('./output/*')

% Import data from corresponding folder
copyfile '../0-Prepare-Data/output/' './input/'

% Tables with menu information and choices
menuTab = readtable('./input/menuTab.csv');
obsTab = readtable('./input/obsTab.csv');

% Set timer
time_start = tic;
save time_start time_start;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Estimation Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Level of background consumption
% recall that payoffs in the menuTab are expressed as dollars
omega = 7.046;

%%% Integration Bounds and number of nodes
rMin = -2  ;
rMax =  2  ;
dMin = -6  ;
dMax =  6  ;
nNodes = 10000;

%%% Initial values
theta_0_pooled = [ ...
    % mu_ra, sig_ra, mu_dr, sig_dr, rho
    -0.052, 0.534, 0.041, 1.72, -0.306
    ];

%%% Parameter bounds
thetaBounds.max.mu_ra  =  1000000;
thetaBounds.min.mu_ra  = -1000000;
thetaBounds.max.sig_ra =  1000000;
thetaBounds.min.sig_ra =  0.00001;
thetaBounds.max.mu_dr  =  1000000;
thetaBounds.min.mu_dr  = -1000000;
thetaBounds.max.sig_dr =  1000000;
thetaBounds.min.sig_dr =  0.00001;
thetaBounds.max.rho    =  0.99999;
thetaBounds.min.rho    = -0.99999;


%%% Settings of gradient-free optimization used to improve initial values
opt_x0 = optimset(...
    'Display', 'iter', ...
    'MaxIter', 50, ...
    'MaxFunEvals', 100000);

%%% Settings of algorithms for optimization used to find MLE
optFminunc_mle = optimoptions('fminunc',...
    'Display','iter-detailed',...
    'OptimalityTolerance',1e-6,...
    'StepTolerance',1e-6,...
    'MaxFunctionEvaluations',5000);

%%% Store in a structure
estInput = struct(...
    'menuTab'          , menuTab           ,...
    'obsTab'           , obsTab            ,...
    'rMin'             , rMin              ,...
    'rMax'             , rMax              ,...
    'dMin'             , dMin              ,...
    'dMax'             , dMax              ,...
    'nNodes'           , nNodes            ,...
    'omega'            , omega             ,...
    'thetaBounds'      , thetaBounds       ,...
    'opt_x0'           , opt_x0            ,...
    'optFminunc_mle'   , optFminunc_mle   );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2) Estimate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial values
estInput.theta_0 = theta_0_pooled;

% MLE estimates pooling all the sample
estOutputPooled = estimate_pooled(estInput);

% Export estimates as CSV table
results_pooled = [
    estOutputPooled.theta_hat(1) ; estOutputPooled.se_theta_hat(1);
    estOutputPooled.theta_hat(2) ; estOutputPooled.se_theta_hat(2);
    estOutputPooled.theta_hat(3) ; estOutputPooled.se_theta_hat(3);
    estOutputPooled.theta_hat(4) ; estOutputPooled.se_theta_hat(4);
    estOutputPooled.theta_hat(5) ; estOutputPooled.se_theta_hat(5);
    estOutputPooled.mom.mle.mu_ra_plus ; estOutputPooled.mom.se.mu_ra_plus;    
    estOutputPooled.mom.mle.median_ra ; estOutputPooled.mom.se.median_ra;
    estOutputPooled.mom.mle.iqr_ra    ; estOutputPooled.mom.se.iqr_ra;
    estOutputPooled.mom.mle.median_dr ; estOutputPooled.mom.se.median_dr;
    estOutputPooled.mom.mle.iqr_dr    ; estOutputPooled.mom.se.iqr_dr;
    estOutputPooled.logLike;
    estOutputPooled.exitFlag;
    estOutputPooled.pdfCheck;
    estOutputPooled.nObs;
    ];

Col_Names = {'RDEU_Pooled'};
Row_Names = {...
    'mu_ra'   , 'se_mu_ra'   , ...
    'sig_ra'  , 'se_sig_ra'  , ...
    'mu_dr'   , 'se_mu_dr'   , ...
    'sig_dr'  , 'se_sig_dr'  , ...
    'rho'     , 'se_rho'     , ...
    'mu_ra_p' , 'se_mu_ra_p' , ...
    'med_ra'  , 'se_med_ra'  , ...
    'iqr_ra'  , 'se_iqr_ra'  , ...
    'med_dr'  , 'se_med_dr'  , ...
    'iqr_dr'  , 'se_iqr_dr'  , ...
    'logLike' , 'exitFlag', 'pdfCheck', 'nObs'};

resTab = table(results_pooled,'RowNames',Row_Names,'VariableNames',Col_Names);
writetable(resTab,'./output/rdeu_pooled_mle_qmc.csv','WriteRowNames',true);

disp(resTab)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3) Report computation time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load time_start time_start;
computation_time = toc(time_start);
toc(time_start);
save('./output/computation_time.mat','computation_time');
delete('./time_start.mat');
