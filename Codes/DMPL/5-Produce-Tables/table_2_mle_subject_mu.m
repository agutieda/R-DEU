% This script:
% - Generate table reporting summary statistics of individual-level estimates
%   of mean risk aversion and discounting using the RDEU and alternative models
%   with DMPL data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0) Initialize and Import Estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;

rdeuTab = readtable('./input/rdeu_subject_mle.csv');
rumLTab = readtable('./input/luce_subject_mle.csv');
rumWTab = readtable('./input/wilcox_subject_mle.csv');
speRTab = readtable('./input/spe_subject_ra.csv','ReadRowNames',true);
speTTab = readtable('./input/spe_subject_dr.csv','ReadRowNames',true);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Produce table of summary statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create matrix of estimates by subject
estMat_mu_ra = [rdeuTab.mu_ra, rumLTab.mu_ra, rumWTab.mu_ra, speRTab.avg_ra];

estMat_mu_dr = [rdeuTab.mu_dr, rumLTab.mu_dr, rumWTab.mu_dr, speTTab.avg_dr];

% Create matrix of summary statistics
sumStats_ra = [
    mean(estMat_mu_ra, 1); ...
    std(estMat_mu_ra, 1); ...
    min(estMat_mu_ra, [], 1); ...
    prctile(estMat_mu_ra, 10, 1); ...
    prctile(estMat_mu_ra, 25, 1); ...
    prctile(estMat_mu_ra, 50, 1); ...
    prctile(estMat_mu_ra, 75, 1); ...
    prctile(estMat_mu_ra, 90, 1); ...
    max(estMat_mu_ra, [], 1);
    corr(estMat_mu_ra, speRTab.avg_ra, "type", "Pearson" )';
    corr(estMat_mu_ra, speRTab.avg_ra, "type", "Kendall" )';
    corr(estMat_mu_ra, speRTab.avg_ra, "type", "Spearman" )'  ];

sumStats_dr = [
    mean(estMat_mu_dr, 1); ...
    std(estMat_mu_dr, 1); ...
    min(estMat_mu_dr, [], 1); ...
    prctile(estMat_mu_dr, 10, 1); ...
    prctile(estMat_mu_dr, 25, 1); ...
    prctile(estMat_mu_dr, 50, 1); ...
    prctile(estMat_mu_dr, 75, 1); ...
    prctile(estMat_mu_dr, 90, 1); ...
    max(estMat_mu_dr, [], 1) ;
    corr(estMat_mu_dr, speTTab.avg_dr, "type", "Pearson" )';
    corr(estMat_mu_dr, speTTab.avg_dr, "type", "Kendall" )';
    corr(estMat_mu_dr, speTTab.avg_dr, "type", "Spearman" )'  ];


% Round everything to nDec decimals
nDec = 3;
sumStats_ra = round(sumStats_ra, nDec);
sumStats_dr = round(sumStats_dr, nDec);

% Define columns
RDEU_ra = sumStats_ra(:,1);
RUML_ra = sumStats_ra(:,2);
RUMW_ra = sumStats_ra(:,3);
SPE_ra  = sumStats_ra(:,4);

RDEU_dr = sumStats_dr(:,1);
RUML_dr = sumStats_dr(:,2);
RUMW_dr = sumStats_dr(:,3);
SPE_dr  = sumStats_dr(:,4);

sep_col = nan(12,1);

% Row labels
rowLabels = {'mean', 'sd', 'min', '10%', '25%', 'median', '75%', ...
    '90%', 'max', 'corr: Pearson', 'corr: Kendall', 'corr: Spearman'};

% Col labels
colLabels = {...
    'RDEU: ra', 'LUCE: ra', 'WILCOX: ra', 'SPE: ra', ...
    '-', ...
    'RDEU: dr', 'LUCE: dr', 'WILCOX: dr', 'SPE: dr'};

% Create table
exportTab = table( ...
    RDEU_ra, RUML_ra, RUMW_ra, SPE_ra, sep_col, ...
    RDEU_dr, RUML_dr, RUMW_dr, SPE_dr, ...
    'RowNames', rowLabels, 'VariableNames', colLabels );

% Eliminate nans
exportTab_CSV = convertvars(exportTab, @isnumeric, @nanblank);

% Export as csv
writetable(exportTab_CSV, './output/table_2.csv', ...
    'WriteRowNames', true);

% Display
disp(exportTab_CSV)