function [logLike, rhoY, rhoY_obs] = loglike_fun(par, estInput)
% Log-likelihood of RUM model with CTB data

% Load
nY = estInput.nY;
nM = estInput.nM;
menuTab = estInput.menuTab;
omega = estInput.omega;

% Auxiliary variables
yList = (0:5:100)';
endowment = 100;

% Recover parameters
mu_ra  = par(1);
mu_dr  = par(2);

% Set tremble probability
nu = eps;

% Check that parameters are in the admissible range
badParameter = ...
    mu_ra  >= estInput.thetaBounds.max.mu_ra  || ...
    mu_ra  <= estInput.thetaBounds.min.mu_ra  || ...
    mu_dr  >= estInput.thetaBounds.max.mu_dr  || ...
    mu_dr  <= estInput.thetaBounds.min.mu_dr  ;

if badParameter == 1

    logLike = nan;

else

    % Vectors describing each menu
    q_t = menuTab.q_t;
    q_tk = menuTab.q_tk;
    t_t = menuTab.t_t;
    t_tk = menuTab.t_tk;
    p_t = menuTab.p_t;
    p_tk = menuTab.p_tk;

    % Initialize
    rhoY = zeros(nM,nY);
    for iMenu = 1:nM

        % DEU of the lottery for each combination of (ra,dr)

        % Payoff of each allocation
        c_t = (endowment-yList)*q_t(iMenu);
        c_tk = yList*q_tk(iMenu);

        % Utility of each allocation
        U_t_1  = ( (c_t  + omega).^(1-mu_ra) - omega.^(1-mu_ra) )./(1-mu_ra);
        U_t_2  = ( (0    + omega).^(1-mu_ra) - omega.^(1-mu_ra) )./(1-mu_ra);
        U_tk_1 = ( (c_tk + omega).^(1-mu_ra) - omega.^(1-mu_ra) )./(1-mu_ra);
        U_tk_2 = ( (0    + omega).^(1-mu_ra) - omega.^(1-mu_ra) )./(1-mu_ra);

        % Expected utility of each allocation
        EU_t = p_t(iMenu) .*U_t_1  +  (1-p_t(iMenu)).*U_t_2;
        EU_tk = p_tk(iMenu).*U_tk_1 + (1-p_tk(iMenu)).*U_tk_2;

        % Discounted expected utility of each allocation
        DEU_t  = exp( -mu_dr.*t_t(iMenu)  ) .* EU_t;
        DEU_tk = exp( -mu_dr.*t_tk(iMenu) ) .* EU_tk;
        DEU  = DEU_t + DEU_tk;

        % RUM probability of choosing each allocation
        rhoY(iMenu,:) = exp(DEU)./sum(exp(DEU));

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

    end

end
