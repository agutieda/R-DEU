% This script:
% - Table 4 reporting ML estimates of the RDEU and alternative models using CB
%   data and pooling all observations in the sample, assuming a baseline level
%   of background consumption of 5 USD

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0) Initialize and Import Estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;

rdeuTab = readtable('./input/rdeu_pooled_mle_baseline.csv','ReadRowNames',true);
rumTab  = readtable('./input/rum_pooled_mle_baseline.csv','ReadRowNames',true);
nlsTab  = readtable('./input/nls_pooled_baseline.csv','ReadRowNames',true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Produce Table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Specify columns
rdeuCol = [...
    rdeuTab{'mu_ra',1}     ; rdeuTab{'se_mu_ra',1}  ;
    nan                    ;
    rdeuTab{'sig_ra',1}    ; rdeuTab{'se_sig_ra',1} ;
    nan                    ;
    rdeuTab{'mu_dr',1}     ; rdeuTab{'se_mu_dr',1}  ;
    nan                    ;
    rdeuTab{'sig_dr',1}    ; rdeuTab{'se_sig_dr',1} ;
    nan                    ;
    rdeuTab{'rho',1}       ; rdeuTab{'se_rho',1}    ;
    nan                    ;    
    rdeuTab{'logLike',1} ] ;

% Specify columns
rumCol = [...
    rumTab{'mu_ra',1}     ; rumTab{'se_mu_ra',1}  ;
    nan                   ;
    nan                   ; nan                   ;
    nan                   ;
    rumTab{'mu_dr',1}     ; rumTab{'se_mu_dr',1}  ;
    nan                   ;
    nan                   ; nan                   ;
    nan                   ;
    nan                   ; nan                   ;
    nan                   ;    
    rumTab{'logLike',1} ] ;

nlsCol = [...
    nlsTab{'mu_ra',1}     ; nlsTab{'se_mu_ra',1}  ;
    nan                   ;
    nan                   ; nan                   ;
    nan                   ;
    nlsTab{'mu_dr',1}     ; nlsTab{'se_mu_dr',1}  ;
    nan                   ;
    nan                   ; nan                   ;
    nan                   ;
    nan                   ; nan                   ;
    nan                   ;    
    nan                ]  ;

% Round
nDec = 3;
rdeuCol = round(rdeuCol,nDec);
rumCol  = round(rumCol,nDec);
nlsCol  = round(nlsCol,nDec);

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
    'rho', 'se_rho',...
    '-----' , ...
    'Log-Like'};

% Col labels
colLab_CSV = {'RDEU','RUM', 'NLS'};

% Create table
exportTab_CSV = table(rdeuCol, rumCol, nlsCol,...
    'VariableNames', colLab_CSV, 'RowNames', rowLab_CSV);

% Eliminate nans
exportTab_CSV = convertvars(exportTab_CSV, @isnumeric, @nanblank);

% Export to CSV
writetable(exportTab_CSV, './output/table_4.csv', ...
    'WriteRowNames', true, 'Delimiter', 'comma', 'FileType', 'text');

% Display table
disp(exportTab_CSV);
