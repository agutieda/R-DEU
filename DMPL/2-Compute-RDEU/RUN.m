% RANDOM DISCOUNTED EXPECTED UTILITY
%
% Estimate R-DEU model using Double Multiple Price Lists (DMPL) data
%
% This script:
% - Estimate risk and time preferences at population and subject level
% - Use Double Multiple Price Lists (DMPL) data from AHLR2008
% - Assume R-DEU model where (ra,dr) follow a bivariate normal distribution
% - Use baseline level of background consumption of 118 DKK
% - Compute standard errors of population estimates clustered at subject level
% - Compute estimates with level of background consumption close to zero
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
rng(1);
addpath([pwd,'/src']);

% Clean folders
delete('./input/*')
% delete('./output/*')

% Import data from corresponding folder
copyfile '../0-Prepare-Data/output/' './input/'

% Import estimates by individual from switching points
copyfile '../1-Compute-SPE/output/spe_subject_ra.csv' './input/spe_subject_ra.csv'
copyfile '../1-Compute-SPE/output/spe_subject_dr.csv' './input/spe_subject_dr.csv'

% Load tables with menu information and choices
menuTab = readtable('./input/menuTab.csv');
obsTab = readtable('./input/obsTab.csv');
riskEst_SP = readtable('./input/spe_subject_ra.csv');
timeEst_SP = readtable('./input/spe_subject_dr.csv');

% Set timer
time_start = tic;
save time_start time_start;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Estimation Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Level of background consumption
% recall that payoffs in the menuTab are expressed as thousands of DKK
omega = 118/1000;

%%% Initial values tested in estimation with pooled data
mu_ra_s0  = mean(riskEst_SP.avg_ra);
mu_dr_s0  = mean(timeEst_SP.avg_dr);
sig_ra_s0 = sqrt(var(riskEst_SP.avg_ra) + mean(riskEst_SP.std_ra.^2) );
sig_dr_s0 = sqrt(var(timeEst_SP.avg_dr) + mean(timeEst_SP.std_dr.^2) );

% Uncomment this line to perform search over a set of initial values
theta_0_pooled = [ ...
    % mu_ra,    sig_ra,    mu_dr,    sig_dr,    rho
    0.50    , 0.50     , 0.20    , 0.20     ,  0; ...
    0.50    , 0.50     , 0.20    , 0.20     ,  0.8; ...
    0.50    , 0.50     , 0.20    , 0.20     , -0.80; ...
    0.00    , 1.00     , 0.10    , 0.10     , -0.80; ...
    -0.50   , 0.50     , 0.10    , 0.10     , -0.80; ...
    0.50    , 0.10     , 0.20    , 0.10     , -0.80; ...
    0.50    , 1.00     , 0.50    , 0.20     ,  0.00; ...
    mu_ra_s0, sig_ra_s0, mu_dr_s0, sig_dr_s0,     0; ...
    ];

%%% Initial values tested for all individuals in estimation with individual data
theta_0_subject = [ ...
    % mu_ra, sig_ra
    0.80  , 0.90 ,  0.10,  0.10 , -0.99 ;
    0.50  , 0.50 ,  0.10,  0.10 , -0.99 ;
    -0.50 , 0.50 ,  0.10,  0.10 , -0.99 ;
    0.00  , 3.00 ,  0.10,  0.10 ,  0.00 ;
    1.50  , 0.10 ,  0.10,  0.10 ,  0.99 ;
    0.00  , 0.10 ,  0.10,  0.10 ,  0.75 ;
    -1.50 , 0.10 ,  0.20 , 0.50 ,  0.50 ;
    1.50  , 1.00 ,  0.20 , 0.80 ,  0.25 ;
    0.00  , 1.00 ,  0.80 , 0.80 ,  0.00 ;
    -1.50 , 1.00 ,  0.20 , 0.20 , -0.25 ;
     0.98 , 1.00 ,  0.20 , 0.20 , -0.50 ;
     0.98 , 0.50 ,  0.20 , 0.20 , -0.75 ;
     0.98 , 1.00 ,  0.10 , 0.10 , -0.99 ;
     0.98 , 0.50 ,  0.10 , 0.10 ,  0.99 ;
     0.80 , 0.10 ,  0.10 , 0.02 ,  0.50 ;
     0.30 , 0.10 ,  0.05 , 0.02 , -0.50 ;
     1.50 , 0.10 ,  0.15 , 0.02 , -0.75 ;
     0.40 , 0.05 ,  0.10 , 0.02 , -0.99 ;
     0.10 , 0.05 ,  0.05 , 0.02 , -0.90 ;
     2.00 , 0.05 ,  0.15 , 0.02 ,  0.00 ;
     0.60 , 0.10 ,  0.05 , 0.01 , -0.90 ;
     2.00 , 0.10 ,  0.05 , 0.01 ,  0.90 ;
     1.50 , 0.60 ,  0.15 , 0.05 , -0.99 ;
     0.40 , 0.60 ,  0.10 , 0.05 ,  0.50 ;
     0.10 , 0.60 ,  0.05 , 0.05 , -0.50 ;
     2.00 , 0.60 ,  0.15 , 0.05 , -0.75 ;
     0.60 , 0.60 ,  0.05 , 0.05 , -0.90 ;
     2.00 , 0.60 ,  0.05 , 0.05 ,  0.50 ;
    ];

% (estimate_subject.m also uses estimates from switching points for each subject 
% as initial values)

%%% Numerical options
nNodes = 401; % Number of nodes used for integration

%%% Parameter bounds
thetaBounds.max.mu_ra  =  1000000  ;
thetaBounds.min.mu_ra  = -1000000  ;
thetaBounds.max.sig_ra =  1000000  ;
thetaBounds.min.sig_ra =  0.00100  ; % Can be lower if we increase nNodes
thetaBounds.max.mu_dr  =  1000000  ;
thetaBounds.min.mu_dr  = -1000000  ;
thetaBounds.max.sig_dr =  1000000  ;
thetaBounds.min.sig_dr =  0.00001  ;
thetaBounds.max.rho    =  0.99999  ;
thetaBounds.min.rho    = -0.99999  ;

%%% Settings of gradient-free optimization used to improve initial values
opt_x0 = optimset(...
    'Display', 'off', ...
    'MaxIter', 2000, ...
    'MaxFunEvals', 100000);

%%% Settings of algorithms for optimization used to find MLE
optFminunc_mle = optimoptions('fminunc',...
    'Display','notify-detailed',...
    'OptimalityTolerance',1e-6,...
    'StepTolerance',1e-6,...
    'MaxFunctionEvaluations',5000);

optFminsearch_mle = optimset(...
    'Display', 'off', ...
    'MaxIter', 10000, ...
    'MaxFunEvals', 100000);

%%% Store in a structure
estInput = struct(...
    'menuTab'          , menuTab           ,...
    'obsTab'           , obsTab            ,...
    'riskEst_SP'       , riskEst_SP        ,...
    'timeEst_SP'       , timeEst_SP        ,...
    'omega'            , omega             ,...
    'nNodes'           , nNodes            ,...    
    'thetaBounds'      , thetaBounds       ,...
    'opt_x0'           , opt_x0            ,...
    'optFminunc_mle'   , optFminunc_mle    ,...
    'optFminsearch_mle' , optFminsearch_mle  );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2) Compute baseline population estimates
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
    estOutputPooled.logLike;
    estOutputPooled.logLike_R; estOutputPooled.logLike_T;
    estOutputPooled.exitFlag; estOutputPooled.pdfCheck;
    estOutputPooled.nObs;
    ];

Col_Names = {'RDEU_Pooled_baseline'};
Row_Names = {...
    'mu_ra' , 'se_mu_ra' , ...
    'sig_ra', 'se_sig_ra', ...
    'mu_dr' , 'se_mu_dr' , ...
    'sig_dr', 'se_sig_dr', ...
    'rho'   , 'se_rho'   , ...
    'logLike' , 'logLike_R' , 'logLike_T',...
    'exitFlag', 'pdfCheck', 'nObs', ...
    };
rTab = table(results_pooled,'RowNames',Row_Names,'VariableNames',Col_Names);
writetable(rTab,'./output/rdeu_pooled_mle.csv','WriteRowNames',true);

% Export predicted choice probabilities by menu
menuID_R = (1:40)';
rho_L1_R = estOutputPooled.rho_L1_R;
rhoObs_L1_R = estOutputPooled.rhoObs_L1_R;
Col_Names = {'menuID', 'rho_L1_R_observed', 'rho_L1_R_predicted'};
rhoRTab = table(menuID_R, rhoObs_L1_R, rho_L1_R, 'VariableNames',Col_Names);
writetable(rhoRTab,'./output/rdeu_pooled_rhoR.csv');

menuID_T = (41:100)';
rho_L1_T = estOutputPooled.rho_L1_T;
rhoObs_L1_T = estOutputPooled.rhoObs_L1_T;
Col_Names = {'menuID', 'rho_L1_T_observed', 'rho_L1_T_predicted'};
rhoTTab = table(menuID_T, rhoObs_L1_T,  rho_L1_T, 'VariableNames',Col_Names);
writetable(rhoTTab,'./output/rdeu_pooled_rhoT.csv');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3) Compute estimates by subject
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial values
estInput.theta_0 = theta_0_subject;

% MLE estimates by subject
estOutputSubject = estimate_subject(estInput);

% Export estimates as CSV table
subjectID    = estOutputSubject.subjectList;
ID           = (1:length(subjectID))';
mu_r_vec     = estOutputSubject.theta_hat(:,1);
sig_r_vec    = estOutputSubject.theta_hat(:,2);
mu_d_vec     = estOutputSubject.theta_hat(:,3);
sig_d_vec    = estOutputSubject.theta_hat(:,4);
rho_vec      = estOutputSubject.theta_hat(:,5);
logLike_vec  = estOutputSubject.logLike;
logLikeR_vec = estOutputSubject.logLike_R;
logLikeT_vec = estOutputSubject.logLike_T;
pdfCheck_vec = estOutputSubject.pdfCheck;
exiFlag_vec  = estOutputSubject.exitFlag;
nObs_vec     = estOutputSubject.nObs;

Col_Names = {'ID', 'subjectID', 'mu_ra', 'sig_ra', 'mu_dr', 'sig_dr', 'rho', ...
    'logLike', 'logLike_R', 'logLike_T', 'exitFlag', 'nObs', 'pdfCheck'};
sTab = table(ID, subjectID, mu_r_vec, sig_r_vec, mu_d_vec, sig_d_vec, rho_vec, ...
    logLike_vec, logLikeR_vec, logLikeT_vec, exiFlag_vec, nObs_vec, pdfCheck_vec, ...
    'VariableNames',Col_Names);
writetable(sTab,'./output/rdeu_subject_mle.csv');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4) Compute population estimates with background consumption close to zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Change level of background consumption
estInput.omega = 1/1000;

% Set initial values
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
    estOutputPooled.logLike;
    estOutputPooled.logLike_R; estOutputPooled.logLike_T;
    estOutputPooled.exitFlag; estOutputPooled.pdfCheck;
    estOutputPooled.nObs;
    ];

Col_Names = {'RDEU_Pooled_w0'};
Row_Names = {...
    'mu_ra' , 'se_mu_ra' , ...
    'sig_ra', 'se_sig_ra', ...
    'mu_dr' , 'se_mu_dr' , ...
    'sig_dr', 'se_sig_dr', ...
    'rho'   , 'se_rho'   , ...
    'logLike' , 'logLike_R' , 'logLike_T',...
    'exitFlag', 'pdfCheck', 'nObs', ...
    };
rTab = table(results_pooled,'RowNames',Row_Names,'VariableNames',Col_Names);
writetable(rTab,'./output/rdeu_pooled_mle_w0.csv','WriteRowNames',true);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5) Report computation time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load time_start time_start;
computation_time = toc(time_start);
toc(time_start);
save('./output/computation_time.mat','computation_time');
delete('./time_start.mat');

% Clean folders
delete('./input/*')