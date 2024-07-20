function [logLike, rhoY, rhoY_obs, pdfCHECK] = loglike_fun(par, estInput)
% Log-likelihood of RDEU model with CTB data

% Load
nY = estInput.nY;
nM = estInput.nM;
chosenA = estInput.chosenA;

% Recover parameters
mu_ra    = par(1); 
sig_ra   = par(2); 
mu_dr    = par(3); 
sig_dr   = par(4); 
mu_pb    = par(5); 
sig_pb   = par(6); 
rho_radr = par(7); 
rho_rapb = par(8); 
rho_drpb = par(9); 

% Set tremble probability
nu = eps;

% Correlation matrix
P = [
        1     ,  rho_radr ,  rho_rapb ;
    rho_radr  ,      1    ,  rho_drpb ;
    rho_rapb  ,  rho_drpb ,     1
    ];

% Check that P is positive semidefinite
eigenvalP = eig(P);
PSDcheck = all(eigenvalP>=0);

% Check that parameters are in the admissible range
badParameter = ...
    mu_ra  >= estInput.thetaBounds.max.mu_ra  || ...
    mu_ra  <= estInput.thetaBounds.min.mu_ra  || ...
    sig_ra >= estInput.thetaBounds.max.sig_ra || ...
    sig_ra <= estInput.thetaBounds.min.sig_ra || ...
    mu_dr  >= estInput.thetaBounds.max.mu_dr  || ...
    mu_dr  <= estInput.thetaBounds.min.mu_dr  || ...
    sig_dr >= estInput.thetaBounds.max.sig_dr || ...
    sig_dr <= estInput.thetaBounds.min.sig_dr || ...
    mu_pb  >= estInput.thetaBounds.max.mu_pb  || ...
    mu_pb  <= estInput.thetaBounds.min.mu_pb  || ...
    sig_pb >= estInput.thetaBounds.max.sig_pb || ...
    sig_pb <= estInput.thetaBounds.min.sig_pb || ...        
    rho_radr >= estInput.thetaBounds.max.rho_radr || ...
    rho_radr <= estInput.thetaBounds.min.rho_radr || ...
    rho_rapb >= estInput.thetaBounds.max.rho_radr || ...
    rho_rapb <= estInput.thetaBounds.min.rho_radr || ...
    rho_drpb >= estInput.thetaBounds.max.rho_drpb || ...
    rho_drpb <= estInput.thetaBounds.min.rho_drpb || ...
    PSDcheck==0;

if badParameter == 1

    logLike = nan;

else

    % Variables used for numerical integration
    rMax = estInput.rMax;
    rMin = estInput.rMin;
    dMax = estInput.dMax;
    dMin = estInput.dMin;
    bMax = estInput.bMax;
    bMin = estInput.bMin;
    rNodes = estInput.rNodes;
    dNodes = estInput.dNodes;
    bNodes = estInput.bNodes;
    intWeights = estInput.intWeights;

    % Compute value of pdf of joint distribution at given points
    PDF_ra = pdf_ra( rNodes , mu_ra , sig_ra , rMin , rMax) ;
    PDF_dr = pdf_dr( dNodes , mu_dr , sig_dr , dMin , dMax) ;
    PDF_pb = pdf_pb( bNodes , mu_pb , sig_pb , bMin , bMax) ;
    CDF_ra = cdf_ra( rNodes , mu_ra , sig_ra , rMin , rMax) ;
    CDF_dr = cdf_dr( dNodes , mu_dr , sig_dr , dMin , dMax) ;
    CDF_pb = cdf_pb( bNodes , mu_pb , sig_pb , bMin , bMax) ;
    copulaC = copulapdf('Gaussian',[CDF_ra,CDF_dr,CDF_pb],P);
    PDF = PDF_ra.*PDF_dr.*PDF_pb.*copulaC;

    % Make sure integration is working
    pdfCHECK = intWeights'*PDF;   

    if abs(pdfCHECK-1)>0.01

        logLike  = nan;
        rhoY     = nan;        
        warning('Integration is not working! Increase number of points!');

    else

        % Create appropiate matrices for integration
        PDF_mat = PDF*ones(1,nY);

        % Compute the probability of allocations 1 to nY for each lottery
        rhoY = nan(nM,nY);
        for iMenu = 1:nM
            integrand_value = chosenA{iMenu}.*PDF_mat;
            rhoY(iMenu,:) =  (intWeights')*integrand_value;
        end

        % Add tremble
        rhoY = rhoY*(1-nu) + nu*(1/nY);

        % Compute the log-likelihood
        idxLogLike  = sub2ind([nM,nY],estInput.obsTab.menuID,estInput.obsTab.Y);
        rhoY_chosen = rhoY(idxLogLike);
        logLike = mean( log(rhoY_chosen) );

    end

    % Additional output
    if nargout > 1

        % Compute observed probability of choosing each menu
        rhoY_obs = nan(nM,nY);
        for jMenu = 1:nM
            for iY = 1:nY
                Y_j = estInput.obsTab.Y(estInput.obsTab.menuID==jMenu);
                rhoY_obs(jMenu,iY) = sum(Y_j==iY)./length(Y_j);
            end
        end

    end

end
