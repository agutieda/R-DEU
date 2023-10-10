% RANDOM DISCOUNTED EXPECTED UTILITY
%
% Estimate NLS model using Convex Budget (CB) data
%
% This script:
% - Estimate risk and time preferences at population level
% - Use CB data from:
%    "Risk Preferences Are Not Time Preferences"
%     by James Andreoni and Charles Sprenger.
%     American Economic Review, (2012)
% - Assume RDEU model where (ra,dr) is fixed and use NLS on FOC
% - Use baseline background consumption of 5 USD
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
omega = 5;

%%% Initial values
theta_0_pooled = [ ...
    % mu_ra, mu_dr
    0.31, 0.26
    ];

%%% Parameter bounds
thetaBounds.max.mu_ra  =  1000000;
thetaBounds.min.mu_ra  = -1000000;
thetaBounds.max.mu_dr  =  1000000;
thetaBounds.min.mu_dr  = -1000000;

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
    'omega'            , omega             ,...
    'thetaBounds'      , thetaBounds       ,...
    'opt_x0'           , opt_x0            ,...
    'optFminunc_mle'   , optFminunc_mle   );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2) Baseline Estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial values
estInput.theta_0 = theta_0_pooled;

% MLE estimates pooling all the sample
estOutputPooled = estimate_pooled(estInput);

% Export estimates as CSV table
results_pooled = [
    estOutputPooled.theta_hat(1) ; estOutputPooled.se_theta_hat(1);
    estOutputPooled.theta_hat(2) ; estOutputPooled.se_theta_hat(2);
    estOutputPooled.exitFlag;
    estOutputPooled.nObs;
    ];

Col_Names = {'NLS_Pooled'};
Row_Names = {...
    'mu_ra'   , 'se_mu_ra'   , ...
    'mu_dr'   , 'se_mu_dr'   , ...
    'exitFlag', 'nObs'};

resTab = table(results_pooled,'RowNames',Row_Names,'VariableNames',Col_Names);
writetable(resTab,'./output/nls_pooled_baseline.csv','WriteRowNames',true);

disp(resTab)


% Compute predicted choice probabilities by menu
a_vec = [0.025:0.05:1, 1];
a_hat = estOutputPooled.a_hat;
a_c = nan(84,1);
rhoY_hat = zeros(84,21);
for jMenu = 1:84

    a_j = a_hat(estInput.obsTab.menuID==jMenu);
    if length(unique(a_j))>1
        warning('Non unique alpha')
    else
        a_j = a_j(1);
        a_c(jMenu) = a_j;
    end

    if a_j<=a_vec(1)
        rhoY_hat(jMenu,1) = 1;
    end

    for iY = 2:21
        if a_j>a_vec(iY-1) &&  a_j<=a_vec(iY)
            rhoY_hat(jMenu,iY) = 1;
        end
    end
end

% Compute observed probability of choosing each menu
rhoY_obs = nan(84,21);
for jMenu = 1:84
    for iY = 1:21
        Y_j = estInput.obsTab.Y(estInput.obsTab.menuID==jMenu);
        rhoY_obs(jMenu,iY) = sum(Y_j==iY)./length(Y_j);
    end
end

menuID = (1:84)';

rhoY_hat_1  = rhoY_hat(:,1 );
rhoY_hat_2  = rhoY_hat(:,2 );
rhoY_hat_3  = rhoY_hat(:,3 );
rhoY_hat_4  = rhoY_hat(:,4 );
rhoY_hat_5  = rhoY_hat(:,5 );
rhoY_hat_6  = rhoY_hat(:,6 );
rhoY_hat_7  = rhoY_hat(:,7 );
rhoY_hat_8  = rhoY_hat(:,8 );
rhoY_hat_9  = rhoY_hat(:,9 );
rhoY_hat_10 = rhoY_hat(:,10);
rhoY_hat_11 = rhoY_hat(:,11);
rhoY_hat_12 = rhoY_hat(:,12);
rhoY_hat_13 = rhoY_hat(:,13);
rhoY_hat_14 = rhoY_hat(:,14);
rhoY_hat_15 = rhoY_hat(:,15);
rhoY_hat_16 = rhoY_hat(:,16);
rhoY_hat_17 = rhoY_hat(:,17);
rhoY_hat_18 = rhoY_hat(:,18);
rhoY_hat_19 = rhoY_hat(:,19);
rhoY_hat_20 = rhoY_hat(:,20);
rhoY_hat_21 = rhoY_hat(:,21);

rhoY_obs_1  = rhoY_obs(:,1 );
rhoY_obs_2  = rhoY_obs(:,2 );
rhoY_obs_3  = rhoY_obs(:,3 );
rhoY_obs_4  = rhoY_obs(:,4 );
rhoY_obs_5  = rhoY_obs(:,5 );
rhoY_obs_6  = rhoY_obs(:,6 );
rhoY_obs_7  = rhoY_obs(:,7 );
rhoY_obs_8  = rhoY_obs(:,8 );
rhoY_obs_9  = rhoY_obs(:,9 );
rhoY_obs_10 = rhoY_obs(:,10);
rhoY_obs_11 = rhoY_obs(:,11);
rhoY_obs_12 = rhoY_obs(:,12);
rhoY_obs_13 = rhoY_obs(:,13);
rhoY_obs_14 = rhoY_obs(:,14);
rhoY_obs_15 = rhoY_obs(:,15);
rhoY_obs_16 = rhoY_obs(:,16);
rhoY_obs_17 = rhoY_obs(:,17);
rhoY_obs_18 = rhoY_obs(:,18);
rhoY_obs_19 = rhoY_obs(:,19);
rhoY_obs_20 = rhoY_obs(:,20);
rhoY_obs_21 = rhoY_obs(:,21);

rhoTab = table( menuID, ...
    rhoY_obs_1 , rhoY_hat_1 , ...
    rhoY_obs_2 , rhoY_hat_2 , ...
    rhoY_obs_3 , rhoY_hat_3 , ...
    rhoY_obs_4 , rhoY_hat_4 , ...
    rhoY_obs_5 , rhoY_hat_5 , ...
    rhoY_obs_6 , rhoY_hat_6 , ...
    rhoY_obs_7 , rhoY_hat_7 , ...
    rhoY_obs_8 , rhoY_hat_8 , ...
    rhoY_obs_9 , rhoY_hat_9 , ...
    rhoY_obs_10, rhoY_hat_10, ...
    rhoY_obs_11, rhoY_hat_11, ...
    rhoY_obs_12, rhoY_hat_12, ...
    rhoY_obs_13, rhoY_hat_13, ...
    rhoY_obs_14, rhoY_hat_14, ...
    rhoY_obs_15, rhoY_hat_15, ...
    rhoY_obs_16, rhoY_hat_16, ...
    rhoY_obs_17, rhoY_hat_17, ...
    rhoY_obs_18, rhoY_hat_18, ...
    rhoY_obs_19, rhoY_hat_19, ...
    rhoY_obs_20, rhoY_hat_20, ...
    rhoY_obs_21, rhoY_hat_21 );

aTab = [menuTab,table(a_c)];

writematrix(rhoY_hat,'./output/nls_rhoY_hat.csv');
writematrix(rhoY_obs,'./output/nls_rhoY_obs.csv');
writetable(rhoTab,'./output/nls_pooled_rhoY_baseline.csv');
writetable(aTab,'./output/nls_pooled_aTab.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3) Report computation time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load time_start time_start;
computation_time = toc(time_start);
toc(time_start);
save('./output/computation_time.mat','computation_time');
delete('./time_start.mat');

% Clean folders
delete('./input/*')