function estInput = compute_thresholds(estInput)
% Function to compute the thresholds that make an individual with RDEU
% preferences indifferent between the two alternatives in a risk menu.
% Also computes the discount rate that, conditional on a value of "r", makes the
% individual indifferent between the two alternatives in a time menu.

menuTab = estInput.menuTab;

% Separate menuTab into risk and time menus
riskMenuTab = menuTab(strcmp(menuTab.taskType, 'Risk'), :);
timeMenuTab = menuTab(strcmp(menuTab.taskType, 'Time'), :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Risk Menus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract info about risk menus and level of background consumption
nRiskMenu = size(riskMenuTab, 1);
omega = estInput.omega;

% Set bounds to find theresholds
if isfield(estInput,'rMin')
    rMin = estInput.rMin;
else
    rMin = -5 ;
end

if isfield(estInput,'rMax')
    rMax = estInput.rMax;
else
    rMax = 5 ;
end

% Loop over risk menus
K_Ar = nan(nRiskMenu, 1);
for iMenu = 1:nRiskMenu
    % Extract menu
    menu = riskMenuTab(iMenu, :);
    % Check there is switching in the interval
    f_l = sign( eu_diff(rMin, menu, omega) );
    f_u = sign( eu_diff(rMax, menu, omega) );
    if f_l == f_u
        K_Ar(iMenu) = nan;
    else
        K_Ar(iMenu) = fzero( @(r) eu_diff(r, menu, omega), [rMin, rMax] );
    end
end

estInput.K_Ar = K_Ar;
estInput.menuTab = menuTab;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time Menus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract info about time menus
nTimeMenu = size(timeMenuTab, 1);
if isfield(estInput,'rNodes')
    nNodes = estInput.nNodes;
    rNodes = estInput.rNodes;
else
    nNodes   = 1 ;
    rNodes   = 0 ;
end

% Loop over risk menus
K_At_r = nan(nTimeMenu, nNodes);

for iMenu = 1:nTimeMenu
    % Extract menu
    menu = timeMenuTab(iMenu, :);
    % Extract info about menu
    x_t  = menu.x1_L1;
    x_tk = menu.x1_L2;
    k    = menu.t_L2 - menu.t_L1;
    num  = (x_tk + omega).^(1 - rNodes) - omega.^(1 - rNodes);
    den  = (x_t  + omega).^(1 - rNodes) - omega.^(1 - rNodes);
    % Compute threshold
    K_At_r(iMenu, :) = (1/k) * log( num./den );
end

estInput.K_At_r = K_At_r;


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Auxiliary Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function util = u_fun(r, x, w)
% Utility function
util = ( (x + w).^(1 - r) - w.^(1 - r)  ) ./ (1 - r);

end


function eu_diff = eu_diff(r, menu, w)
% Function that computes difference in EU between lotteries

% Utility in each outcome
U1_L1 = u_fun(r, menu.x1_L1, w);
U2_L1 = u_fun(r, menu.x2_L1, w);
U1_L2 = u_fun(r, menu.x1_L2, w);
U2_L2 = u_fun(r, menu.x2_L2, w);

% Expected utility of each lottery
EU_L1 = menu.p_L1 .* U1_L1 + (1 - menu.p_L1) .* U2_L1;
EU_L2 = menu.p_L2 .* U1_L2 + (1 - menu.p_L2) .* U2_L2;

% Difference in EU
eu_diff = EU_L1 - EU_L2;

end
