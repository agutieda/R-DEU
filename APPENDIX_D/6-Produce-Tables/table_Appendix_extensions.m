% RANDOM DISCOUNTED EXPECTED UTILITY
%
% Generate table reported in the supplementary material
%
% This script:
% - Generate Appendix Table reporting population-level estimates of risk and time
%   preferences using the R-DEU model with exponential and hyperbolic discounting
%   under different assumptions about the distribution of preferences.
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
copyfile '../1-Compute-Baseline/output/' './input/'
copyfile '../2-Compute-QMC/output/' './input/'
copyfile '../3-Compute-QMC-Constrained/output/' './input/'
copyfile '../4-Compute-HD/output/' './input/'
copyfile '../5-Compute-HD-Constrained/output/' './input/'

% Load estimated parameters under each specification
baseTab = readtable('./input/rdeu_pooled_mle_baseline.csv','ReadRowNames',true);
qmcTab  = readtable('./input/rdeu_pooled_mle_qmc.csv','ReadRowNames',true);
qmcXTab = readtable('./input/rdeu_pooled_mle_qmc_constrained.csv','ReadRowNames',true);
hbdTab  = readtable('./input/rdeu_pooled_mle_hd.csv','ReadRowNames',true);
hbdXTab = readtable('./input/rdeu_pooled_mle_hd_constrained.csv','ReadRowNames',true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Produce Table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify columns
baseCol = [...
    baseTab{'mu_ra',1}     ; baseTab{'se_mu_ra',1}  ;
    nan                    ;
    baseTab{'iqr_ra',1}    ; baseTab{'se_iqr_ra',1} ;
    nan                    ;
    baseTab{'mu_dr',1}     ; baseTab{'se_mu_dr',1}  ;
    nan                    ;
    baseTab{'iqr_dr',1}    ; baseTab{'se_iqr_dr',1} ;
    nan                    ;
    nan                    ;
    nan                    ; nan                    ;
    nan                    ;
    nan                    ; nan                    ;
    baseTab{'rho',1}       ; baseTab{'se_rho',1}    ;
    nan                    ;
    nan                    ; nan                    ;
    nan                    ;
    nan                    ; nan                    ;
    nan                    ;
    baseTab{'logLike',1} ] ;

qmcCol = [...
    qmcTab{'med_ra',1}     ; qmcTab{'se_med_ra',1}  ;
    nan                    ;
    qmcTab{'iqr_ra',1}     ; qmcTab{'se_iqr_ra',1}  ;
    nan                    ;
    qmcTab{'med_dr',1}     ; qmcTab{'se_med_dr',1}  ;
    nan                    ;
    qmcTab{'iqr_dr',1}     ; qmcTab{'se_iqr_dr',1}  ;
    nan                    ;
    nan                    ; nan                    ;
    nan                    ;
    nan                    ; nan                    ;
    nan                    ;
    qmcTab{'rho',1}        ; qmcTab{'se_rho',1}     ;
    nan                    ;
    nan                    ; nan                    ;
    nan                    ;
    nan                    ; nan                    ;
    nan                    ;
    qmcTab{'logLike',1} ] ;

qmcXCol = [...
    qmcXTab{'med_ra',1}    ; qmcXTab{'se_med_ra',1} ;
    nan                    ;
    qmcXTab{'iqr_ra',1}    ; qmcXTab{'se_iqr_ra',1} ;
    nan                    ;
    qmcXTab{'med_dr',1}    ; qmcXTab{'se_med_dr',1} ;
    nan                    ;
    qmcXTab{'iqr_dr',1}    ; qmcXTab{'se_iqr_dr',1} ;
    nan                    ;
    nan                    ; nan                    ;
    nan                    ;
    nan                    ; nan                    ;
    nan                    ;
    qmcXTab{'rho',1}       ; qmcXTab{'se_rho',1}    ;
    nan                    ;
    nan                    ; nan                    ;
    nan                    ;
    nan                    ; nan                    ;
    nan                    ;
    qmcXTab{'logLike',1} ] ;

hbdCol = [...
    hbdTab{'med_ra',1}     ; hbdTab{'se_med_ra',1}  ;
    nan                    ;
    hbdTab{'iqr_ra',1}     ; hbdTab{'se_iqr_ra',1}  ;
    nan                    ;
    hbdTab{'med_dr',1}     ; hbdTab{'se_med_dr',1}  ;
    nan                    ;
    hbdTab{'iqr_dr',1}     ; hbdTab{'se_iqr_dr',1}  ;
    nan                    ;
    hbdTab{'med_pb',1}     ; hbdTab{'se_med_pb',1}  ;
    nan                    ;
    hbdTab{'iqr_pb',1}     ; hbdTab{'se_iqr_pb',1}  ;
    nan                    ;
    hbdTab{'rho_radr',1}   ; hbdTab{'se_rho_radr',1};
    nan                    ;
    hbdTab{'rho_rapb',1}   ; hbdTab{'se_rho_rapb',1};
    nan                    ;
    hbdTab{'rho_drpb',1}   ; hbdTab{'se_rho_drpb',1};
    nan                    ;
    hbdTab{'logLike',1} ]  ;

hbdXCol = [...
    hbdXTab{'med_ra',1}     ; hbdXTab{'se_med_ra',1}  ;
    nan                    ;
    hbdXTab{'iqr_ra',1}     ; hbdXTab{'se_iqr_ra',1}  ;
    nan                    ;
    hbdXTab{'med_dr',1}     ; hbdXTab{'se_med_dr',1}  ;
    nan                    ;
    hbdXTab{'iqr_dr',1}     ; hbdXTab{'se_iqr_dr',1}  ;
    nan                    ;
    hbdXTab{'med_pb',1}     ; hbdXTab{'se_med_pb',1}  ;
    nan                    ;
    hbdXTab{'iqr_pb',1}     ; hbdXTab{'se_iqr_pb',1}  ;
    nan                    ;
    hbdXTab{'rho_radr',1}   ; hbdXTab{'se_rho_radr',1};
    nan                    ;
    hbdXTab{'rho_rapb',1}   ; hbdXTab{'se_rho_rapb',1};
    nan                    ;
    hbdXTab{'rho_drpb',1}   ; hbdXTab{'se_rho_drpb',1};
    nan                    ;
    hbdXTab{'logLike',1} ]  ;

% Round
nDec = 4;
baseCol = round(baseCol,nDec);
qmcCol  = round(qmcCol,nDec);
qmcXCol = round(qmcXCol,nDec);
hbdCol  = round(hbdCol,nDec);
hbdXCol = round(hbdXCol,nDec);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2) Export to CSV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Column labels
rowLab_CSV = {...
    'mu_ra', 'se_mu_ra', ...
    '-' , ...
    'sigma_ra', 'se_sigma_ra', ...
    '--' , ...
    'mu_dr', 'se_mu_dr', ...
    '---' , ...
    'sigma_dr', 'se_sigma_dr',...
    '----' , ...
    'mu_pb', 'se_mu_pb', ...
    '-----' , ...
    'sigma_pb', 'se_sigma_pb',...
    '------' , ...
    'rho_radr', 'se_rho_radr',...
    '-------' , ...
    'rho_rapb', 'se_rho_rapb',...
    '--------' , ...
    'rho_drpb', 'se_rho_drpb',...
    '---------' , ...
    'Log-Like'};

% Col labels
colLab_CSV = {...
    'Baseline', 'QMC', 'QMC-Constrained', 'Hyperbolic-Discounting' 'Hyperbolic-Discounting-Constrained'};

% Create table
exportTab_CSV = table( ...
    baseCol, qmcCol, qmcXCol, hbdCol, hbdXCol, ...
    'VariableNames', colLab_CSV, 'RowNames', rowLab_CSV);

% Eliminate nans
exportTab_CSV = convertvars(exportTab_CSV, @isnumeric, @nanblank);

% Export to CSV
writetable(exportTab_CSV, './output/table_appendix_extensions.csv', ...
    'WriteRowNames', true, 'Delimiter', 'comma', 'FileType', 'text');

% Display table
disp(exportTab_CSV);

% Clean folders
delete('./input/*')
