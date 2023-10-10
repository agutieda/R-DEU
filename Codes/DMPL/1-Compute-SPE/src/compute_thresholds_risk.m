function K_Ar = compute_thresholds_risk(menuTab, omega)
% Function to compute the thresholds that make an individual with RDEU
% preferences indifferent between the two alternatives in a risk menu.
% Also computes the discount rate that, conditional on a value of "r", makes the
% individual indifferent between the two alternatives in a time menu.

% Select risk menus in menuTab
riskMenuTab = menuTab(strcmp(menuTab.taskType, 'Risk'), :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Risk Menus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract info about risk menus and level of background consumption
nRiskMenu = size(riskMenuTab, 1);

% Set bounds to find theresholds
rMin = -10;
rMax =  10;

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
