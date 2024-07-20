function [SSE,a_hat] = sse_fun(par, estInput)
% Objective function of model with CTB data to estimate using NLS

% Load
menuTab = estInput.menuTab;
omega = estInput.omega;

% Auxiliary variables
endowment = 100;

% Recover parameters
mu_ra  = par(1);
mu_dr  = par(2);

% Check that parameters are in the admissible range
badParameter = ...
    mu_ra  >= estInput.thetaBounds.max.mu_ra  || ...
    mu_ra  <= estInput.thetaBounds.min.mu_ra  || ...
    mu_dr  >= estInput.thetaBounds.max.mu_dr  || ...
    mu_dr  <= estInput.thetaBounds.min.mu_dr  ;

if badParameter == 1

    SSE = nan;

else

    % Vectors describing each menu
    q_t = menuTab.q_t;
    q_tk = menuTab.q_tk;
    p_t = menuTab.p_t;
    p_tk = menuTab.p_tk;
    k = menuTab.k;

    % auxiliary variables
    p_ratio = p_tk./p_t;
    R_k = q_tk./q_t;
    m = endowment*q_tk;

    % Consumption predicted for each menu
    auxC = ( exp(-mu_dr.*k).*p_ratio.*R_k ).^(-1/mu_ra);
    beta1 = auxC./(1+auxC.*R_k);
    beta2 = 1./(1+auxC.*R_k);
    c_t_Menu = beta1.*(m + omega) - beta2.*omega;

    % Consumption predicted for each observation
    a_hat = max(min( 1 - c_t_Menu./(endowment.*q_t), 1 ),0) ;
    a_hat = a_hat(estInput.obsTab.menuID);
    c_t_hat = c_t_Menu(estInput.obsTab.menuID);

    % Consumption observed in each menu
    c_t_obs = estInput.obsTab.c_t;

    % Sum of squared errors
    SSE = mean( (c_t_obs - c_t_hat).^2 );    

end

end
