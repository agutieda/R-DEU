function K_At_r = compute_thresholds_time(menuTab, omega, varargin)
% Function to compute the thresholds that make an individual with RDEU
% preferences indifferent between the two alternatives in a risk menu.
% Also computes the discount rate that, conditional on a value of "r", makes the
% individual indifferent between the two alternatives in a time menu.

% Select time menus in menuTab
timeMenuTab = menuTab(strcmp(menuTab.taskType, 'Time'), :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time Menus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract info about time menus
nTimeMenu = size(timeMenuTab, 1);
if nargin>2
    rac = varargin{1};
    nNodes = length(rac);
else
    rac = 0;
    nNodes = 1;
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
    num  = (x_tk + omega).^(1 - rac) - omega.^(1 - rac);
    den  = (x_t  + omega).^(1 - rac) - omega.^(1 - rac);
    % Compute threshold
    K_At_r(iMenu, :) = (1/k) * log( num./den );
end


end
