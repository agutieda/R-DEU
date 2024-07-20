% This script:
% - Generate table reporting population-level estimates of risk aversion and
%   discounting under the RDEU and alternative models using DMPL data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0) Initialize and Import Estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;

rdeuTab = readtable('./input/rdeu_pooled_mle.csv','ReadRowNames',true);
rumLTab = readtable('./input/luce_pooled_mle.csv','ReadRowNames',true);
rumWTab = readtable('./input/wilcox_pooled_mle.csv','ReadRowNames',true);
spRTab = readtable('./input/spe_pooled_ra.csv');
spTTab = readtable('./input/spe_pooled_dr.csv');

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
    rdeuTab{'logLike_R',1} ; rdeuTab{'logLike_T',1} ;
    rdeuTab{'logLike',1} ] ;

rumLCol = [...
    rumLTab{'mu_ra',1}     ; rumLTab{'se_mu_ra',1}   ;
    nan                    ;
    rumLTab{'sig_ra',1}    ; rumLTab{'se_sig_ra',1}  ;
    nan                    ;
    rumLTab{'mu_dr',1}     ; rumLTab{'se_mu_dr',1}   ;
    nan                    ;
    rumLTab{'sig_dr',1}    ; rumLTab{'se_sig_dr',1}  ;
    nan                    ;
    nan                    ; nan                     ;
    nan                    ;
    rumLTab{'logLike_R',1} ;  rumLTab{'logLike_T',1} ;
    rumLTab{'logLike',1} ] ;

rumWCol = [...
    rumWTab{'mu_ra',1}     ; rumWTab{'se_mu_ra',1}   ;
    nan                    ;
    rumWTab{'sig_ra',1}    ; rumWTab{'se_sig_ra',1}  ;
    nan                    ;
    rumWTab{'mu_dr',1}     ; rumWTab{'se_mu_dr',1}   ;
    nan                    ;
    rumWTab{'sig_dr',1}    ; rumWTab{'se_sig_dr',1}  ;
    nan                    ;
    nan                    ; nan                     ;
    nan                    ;
    rumWTab{'logLike_R',1} ;  rumWTab{'logLike_T',1} ;
    rumWTab{'logLike',1} ] ;

speCol = [...
    spRTab.mu_ra  ; nan  ;
    nan           ;
    spRTab.sig_ra ; nan  ;
    nan           ;
    spTTab.mu_dr  ; nan  ;
    nan           ;
    spTTab.sig_dr ; nan  ;
    nan           ;
    spTTab.rho    ; nan  ;
    nan           ;
    nan           ; nan  ;
    nan   ];

% Round
nDec = 3;
rdeuCol  = round(rdeuCol,nDec);
rumLCol  = round(rumLCol,nDec);
rumWCol  = round(rumWCol,nDec);
speCol    = round(speCol,nDec);

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
    'Log-Like (Risk Menus)', ...
    'Log-Like (Time Menus)', ...
    'Log-Like'};

% Col labels
colLab_CSV = {...
    'RDEU', 'LUCE', 'WILCOX', 'SPE' };

% Create table
exportTab_CSV = table( ...
    rdeuCol, rumLCol, rumWCol, speCol, ...
    'VariableNames', colLab_CSV, 'RowNames', rowLab_CSV);

% Eliminate nans
exportTab_CSV = convertvars(exportTab_CSV, @isnumeric, @nanblank);

% Export to CSV
writetable(exportTab_CSV, './output/table_1.csv', ...
    'WriteRowNames', true, 'Delimiter', 'comma', 'FileType', 'text');

% Display table
disp(exportTab_CSV);
