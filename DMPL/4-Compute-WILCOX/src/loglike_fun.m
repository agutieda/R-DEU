function [logLike, logLike_R, logLike_T, rho_L1_R, rho_L1_T, ...
    rhoObs_L1_R, rhoObs_L1_T] = loglike_fun(theta_hat, estInput)
% Log-likelihood function of RUM model with DMPL data
% Uses the same Luce specification as in AHLR 2008

% Unpack parameters
mu_ra  = theta_hat(1); % risk aversion
sig_ra = theta_hat(2); % SD of the utility shocks in risk menus
mu_dr  = theta_hat(3); % discount rate
sig_dr = theta_hat(4); % of the utility shocks in time menus

% Check that parameters are in the admissible range
badParameter = ...
    mu_ra  >= estInput.thetaBounds.max.mu_ra  || ...
    mu_ra  <= estInput.thetaBounds.min.mu_ra  || ...
    sig_ra >= estInput.thetaBounds.max.sig_ra || ...
    sig_ra <= estInput.thetaBounds.min.sig_ra || ...
    mu_dr  >= estInput.thetaBounds.max.mu_dr  || ...
    mu_dr  <= estInput.thetaBounds.min.mu_dr  || ...
    sig_dr >= estInput.thetaBounds.max.sig_dr || ...
    sig_dr <= estInput.thetaBounds.min.sig_dr ;

% Check coverage of integral bounds
if badParameter == 1

    logLike = nan;

else

    obsTab = estInput.obsTab;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute prob. of choosing safe alternative (L1) in each risk menu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    menuR = estInput.menuTab(1:40,:);
    omega = estInput.omega;

    % Utility in each outcome
    U1_L1 = ( (omega + menuR.x1_L1).^(1 - mu_ra)  )./(1-mu_ra);
    U2_L1 = ( (omega + menuR.x2_L1).^(1 - mu_ra)  )./(1-mu_ra);
    U1_L2 = ( (omega + menuR.x1_L2).^(1 - mu_ra)  )./(1-mu_ra);
    U2_L2 = ( (omega + menuR.x2_L2).^(1 - mu_ra)  )./(1-mu_ra);

    % Expected utility of each lottery
    DEU_L1 = menuR.p_L1.*U1_L1 + (1 - menuR.p_L1).*U2_L1;
    DEU_L2 = menuR.p_L2.*U1_L2 + (1 - menuR.p_L2).*U2_L2;

    % Maximum and minimum possible utility
    uMax  = ( (omega + estInput.cMax).^(1 - mu_ra)  )./(1-mu_ra);
    uMin  = ( (omega + estInput.cMin).^(1 - mu_ra)  )./(1-mu_ra);

    % Normalized difference of expected utilities
    DEU_DIFF = ((DEU_L1 - DEU_L2)./abs(uMax - uMin))./sig_ra;

    % Compute the probability of choosing L1 in each risky menu
    rho_L1_R = normcdf(DEU_DIFF);

    % Add tremble probability to handle cases where rho==1 or rho==0
    nu = 1e-12;
    rho_L1_R = rho_L1_R*(1 - nu) + 0.5*nu;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Compute prob. of choosing early payoff (L1) in each time menu
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    menuT = estInput.menuTab(41:100,:);

    % Utility in each outcome
    U_L1 = ( (omega + menuT.x1_L1).^(1 - mu_ra)  )./(1-mu_ra);
    U_L2 = ( (omega + menuT.x1_L2).^(1 - mu_ra)  )./(1-mu_ra);
    U_w =  ( omega.^(1 - mu_ra)  )./(1-mu_ra);

    % Discounted expected utility
    DEU_L1 = exp(-mu_dr*menuT.t_L1).*U_L1 + exp(-mu_dr*menuT.t_L2).*U_w;
    DEU_L2 = exp(-mu_dr*menuT.t_L1).*U_w  + exp(-mu_dr*menuT.t_L2).*U_L2;

    % Normalized difference of expected utilities
    DEU_DIFF = (DEU_L1 - DEU_L2)./sig_dr;

    % Compute the probability of choosing L1 in each risky menu
    rho_L1_T = normcdf(DEU_DIFF);

    % Add tremble probability to handle cases where rho==1 or rho==0
    nu = 1e-12;
    rho_L1_T = rho_L1_T*(1 - nu) + 0.5*nu;

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

     % Check log-likelihood is real
     if ~isreal(logLike)
        logLike = nan;
    end

    % Additional output
    if nargout > 1

        % Compute observed probability of choosing each risk menu
        menuList_R = 1:40;
        nMenus_R = length(menuList_R);
        rhoObs_L1_R = nan(nMenus_R,1);
        for jMenu = 1:nMenus_R
            Y_j = obsTabR{obsTabR.menuID==menuList_R(jMenu),{'Y'}};
            nObs_j = length(Y_j);
            rhoObs_L1_R(jMenu) = sum((Y_j == 1) + 0.5*(Y_j == 0))/nObs_j;
        end

        % Compute observed probability of choosing each time menu
        menuList_T = 41:100;
        nMenus_T = length(menuList_T);
        rhoObs_L1_T = nan(nMenus_T,1);
        for jMenu = 1:nMenus_T
            Y_j = obsTabT{obsTabT.menuID==menuList_T(jMenu),{'Y'}};
            nObs_j = length(Y_j);
            rhoObs_L1_T(jMenu) = sum((Y_j == 1) + 0.5*(Y_j == 0))/nObs_j;
        end

    end

end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = invlogit(x)
% Inverse logit function
p = exp(x)./(1 + exp(x));
end