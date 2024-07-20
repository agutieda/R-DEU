% This script:
% - Generate first two columns of Appendix Table reporting population-level 
%   estimates of risk and time preferences uisng the RDEU model under different 
%   levels of background consumption.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0) Initialize and Import Estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;

rdeu_Tab = readtable('./input/rdeu_pooled_mle.csv', 'ReadRowNames',true);
rdeu_w0_Tab = readtable('./input/rdeu_pooled_mle_w0.csv', 'ReadRowNames',true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Produce Table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify columns
rdeuCol = [...
    rdeu_Tab{'mu_ra',1}     ; rdeu_Tab{'se_mu_ra',1}  ;
    nan                     ;
    rdeu_Tab{'sig_ra',1}    ; rdeu_Tab{'se_sig_ra',1} ;
    nan                     ;
    rdeu_Tab{'mu_dr',1}     ; rdeu_Tab{'se_mu_dr',1}  ;
    nan                     ;
    rdeu_Tab{'sig_dr',1}    ; rdeu_Tab{'se_sig_dr',1} ;
    nan                     ;
    rdeu_Tab{'rho',1}       ; rdeu_Tab{'se_rho',1}   ];

% Specify columns
rdeu_w0_Col = [...
    rdeu_w0_Tab{'mu_ra',1}     ; rdeu_w0_Tab{'se_mu_ra',1}  ;
    nan                    ;
    rdeu_w0_Tab{'sig_ra',1}    ; rdeu_w0_Tab{'se_sig_ra',1} ;
    nan                    ;
    rdeu_w0_Tab{'mu_dr',1}     ; rdeu_w0_Tab{'se_mu_dr',1}  ;
    nan                    ;
    rdeu_w0_Tab{'sig_dr',1}    ; rdeu_w0_Tab{'se_sig_dr',1} ;
    nan                    ;
    rdeu_w0_Tab{'rho',1}       ; rdeu_w0_Tab{'se_rho',1} ] ;

% Round
nDec = 3;
rdeuCol = round(rdeuCol,nDec);
rdeu_w0_Col = round(rdeu_w0_Col,nDec);

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
    'rho', 'se_rho'};

% Col labels
colLab_CSV = {...
    'Baseline', 'omega = 0'};

% Create table
exportTab_CSV = table( ...
    rdeuCol, rdeu_w0_Col, ...
    'VariableNames', colLab_CSV, 'RowNames', rowLab_CSV);

% Eliminate nans
exportTab_CSV = convertvars(exportTab_CSV, @isnumeric, @nanblank);

% Export to CSV
writetable(exportTab_CSV, './output/table_appendix_omega_c12.csv', ...
    'WriteRowNames', true, 'Delimiter', 'comma', 'FileType', 'text');

% Display table
disp(exportTab_CSV);
