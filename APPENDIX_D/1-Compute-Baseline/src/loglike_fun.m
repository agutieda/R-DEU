function [logLike, rhoY, rhoY_obs, pdfCheck] = loglike_fun(par, estInput)
% Log-likelihood of RDEU model with CTB data

% Load
nY = estInput.nY;
nM = estInput.nM;
nA = estInput.nA;

% Recover parameters
mu_ra  = par(1);
sig_ra = par(2);
mu_dr  = par(3);
sig_dr = par(4);
rho    = par(5);

% Set tremble probability
nu = eps;

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
    rho    >= estInput.thetaBounds.max.rho    || ...
    rho    <= estInput.thetaBounds.min.rho;

if badParameter == 1

    logLike = nan;

else

    % Variables used for numerical integration
    K_convex = estInput.K_convex;
    K_concave = estInput.K_concave;
    ra_convex = estInput.ra_convex;
    ra_concave = estInput.ra_concave;
    intWeights_concave = estInput.intWeights_concave;
    intWeights_convex = estInput.intWeights_convex;

    % Integration weights and conditional moments for convex and concave part
    w_phi_ra_cvx = normpdf(ra_convex, mu_ra, sig_ra).*intWeights_convex;
    w_phi_ra_ccv = normpdf(ra_concave, mu_ra, sig_ra).*intWeights_concave;
    mu_dr_ra_cvx = mu_dr + rho*sig_dr*(ra_convex-mu_ra)./sig_ra;
    mu_dr_ra_ccv = mu_dr + rho*sig_dr*(ra_concave-mu_ra)./sig_ra;
    sig_dr_ra = sig_dr*sqrt(1-rho^2);

    % Compute the probability of observing each allocation for each menu
    rhoY = zeros(nM,nY);
    for iMenu = 1:nM

        %%% Concave utilities
        rhoY_CCV = zeros(1,nA);
        P_am = 0;
        for iA = 1:nA
            % Phi_dr = 1 - normcdf(K_concave(:,iMenu,iA), mu_dr_ra_ccv, sig_dr_ra);
            z_r = ( K_concave(:,iMenu,iA) - mu_dr_ra_ccv)./sig_dr_ra;
            Phi_dr= 1 -  0.5*erfc(-z_r./sqrt(2)); % Faster than using normcdf
            P_a = w_phi_ra_ccv'*Phi_dr;
            rhoY_CCV(iA) = P_a - P_am;
            P_am = P_a;
        end

        %%% Convex utilities
        % Phi_dr = 1 - normcdf(K_convex(:,iMenu), mu_dr_ra_cvx, sig_dr_ra);
        z_x = (K_convex(:,iMenu) - mu_dr_ra_cvx)./sig_dr_ra;
        Phi_dr= 1 -  0.5*erfc(-z_x./sqrt(2)); % Faster than using normcdf
        rhoY_CVX = w_phi_ra_cvx'*Phi_dr;

        % Probability of allocating all tokens to early date (alpha=0)
        rhoY(iMenu,1) = rhoY_CVX(1) + rhoY_CCV(1);

        % Probability of allocating tokens in between
        rhoY(iMenu,2:(nY-1)) = rhoY_CCV(2:end);

        % Probability of allocating all tokens to delayed date (alpha=1)
        rhoY(iMenu,nY) = 1 - sum(rhoY(iMenu,:));


    end

    % Add tremble
    rhoY = rhoY*(1-nu) + nu*(1/nY);

    % Compute the log-likelihood
    idxLogLike  = sub2ind([nM,nY],estInput.obsTab.menuID,estInput.obsTab.Y);
    rhoY_chosen = rhoY(idxLogLike);
    logLike = mean( log(rhoY_chosen) );

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

        % Check
        pdfCheck = sum(w_phi_ra_cvx) + sum(w_phi_ra_ccv);

    end

end
