function [logLike, logLike_R, logLike_T, pdfCheck, rho_L1_R, rho_L1_T, ...
    rhoObs_L1_R, rhoObs_L1_T] = loglike_fun(theta_hat, estInput)
% Log-likelihood function of RDEU model with DMPL data

% Unpack parameters
mu_ra  = theta_hat(1); % Mean of the distribution of risk aversion
sig_ra = theta_hat(2); % SD of the distribution of risk aversion
mu_dr  = theta_hat(3); % Mean of the distribution of discount rate
sig_dr = theta_hat(4); % SD of the distribution of discount rate

% Correlation coefficient btw ra and dr
if length(theta_hat) == 5
    rho = theta_hat(5);
else
    rho = 0;
end

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

% Check coverage of integral bounds
if badParameter == 1

    logLike = nan;

else

    obsTab = estInput.obsTab;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute prob. of choosing safe alternative (L1) in each risk menu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % rho_L1_R = 1 - normcdf(estInput.K_Ar, mu_ra, sig_ra);
    z_r = ( estInput.K_Ar - mu_ra)./sig_ra;
    rho_L1_R = 1 - 0.5*erfc(-z_r./sqrt(2));

    % Add tremble probability to handle cases where rho==1 or rho ==0
    nu = 1e-12;
    rho_L1_R = rho_L1_R*(1 - nu) + 0.5*nu;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute prob. of choosing early payoff (L1) in each time menu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Integration weights
    w_phi_r = normpdf(estInput.rNodes, mu_ra, sig_ra).*estInput.rWeights;

    % Conditional mean and sd of dr given ra
    mu_dr_cond = mu_dr + sig_dr * rho * (estInput.rNodes - mu_ra)./sig_ra;
    sig_dr_cond = sig_dr * sqrt(1 - rho^2);

    % Probability of choosing the early payoff in each menu
    nTimeMenus = size(estInput.K_At_r,1);
    rho_L1_T = nan(nTimeMenus,1);
    for iMenu = 1:nTimeMenus

        % Phi_dr_cond = normcdf(K_At_r(iMenu,:)', mu_dr_cond, sig_dr_cond);

        % Faster than using normcdf
        z_r = ( estInput.K_At_r(iMenu,:)' - mu_dr_cond)./sig_dr_cond;
        Phi_dr_cond = 0.5*erfc(-z_r./sqrt(2));
        rho_L1_T(iMenu) = 1 - Phi_dr_cond'*w_phi_r;

    end

    % Add tremble probability to handle cases where Y=1
    nu = 1e-12;
    rho_L1_T = rho_L1_T * (1 - nu) + 0.5 * nu;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute log-likehood
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Separate obsTab into risk and time menus
    obsTabR = obsTab(strcmp(obsTab.taskType, 'Risk'), :);
    obsTabT = obsTab(strcmp(obsTab.taskType, 'Time'), :);

    % Log-likelihood of risk menus
    P_L1_R = rho_L1_R(obsTabR.menuID);
    logLike_L1_R = log(P_L1_R).*( (obsTabR.Y == 1) + 0.5*(obsTabR.Y == 0) );
    logLike_L2_R = log(1 - P_L1_R).*( (obsTabR.Y == 2) + 0.5*(obsTabR.Y == 0) );
    logLike_R = mean(logLike_L1_R + logLike_L2_R);

    % Log-likelihood of time menus
    P_L1_T = rho_L1_T(obsTabT.menuID - 40);
    logLike_L1_T = log(P_L1_T).*( (obsTabT.Y == 1) + 0.5*(obsTabT.Y == 0) );
    logLike_L2_T = log(1 - P_L1_T).*( (obsTabT.Y == 2) + 0.5*(obsTabT.Y == 0) );
    logLike_T = mean(logLike_L1_T + logLike_L2_T);

    % Log-likelihood of the whole dataset
    nR = length(P_L1_R);
    nT = length(P_L1_T);
    weightR = nR/(nR+nT);
    logLike = weightR*logLike_R + (1-weightR)*logLike_T;

    % Check integration bounds
    pdfCheck = sum(w_phi_r);

    % Check log-likelihood is real
    if ~isreal(logLike)
        logLike = nan;
    end

    % Additional output
    if nargout > 1

        % Compute observed probability of choosing each risk menu
        nMenus_R = size(estInput.K_Ar,1);
        menuList_R = 1:40;
        rhoObs_L1_R = nan(nMenus_R,1);
        for jMenu = 1:nMenus_R
            Y_j = obsTabR{obsTabR.menuID==menuList_R(jMenu),{'Y'}};
            nObs_j = length(Y_j);
            rhoObs_L1_R(jMenu) = sum((Y_j == 1) + 0.5*(Y_j == 0))/nObs_j;
        end

        % Compute observed probability of choosing each time menu
        nMenus_T = size(estInput.K_At_r,1);
        menuList_T = 41:100;
        rhoObs_L1_T = nan(nMenus_T,1);
        for jMenu = 1:nMenus_T
            Y_j = obsTabT{obsTabT.menuID==menuList_T(jMenu),{'Y'}};
            nObs_j = length(Y_j);
            rhoObs_L1_T(jMenu) = sum((Y_j == 1) + 0.5*(Y_j == 0))/nObs_j;
        end

    end

end


end