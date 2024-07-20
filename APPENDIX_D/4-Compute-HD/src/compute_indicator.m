function estInput = compute_indicator(estInput)
% Function to compute indicator function taking values of one when an alternative
% maximizes DEU in a menu, and zero otherwise

% Menus
q_t  = estInput.menuTab.q_t  ;
q_tk = estInput.menuTab.q_tk ;
p_t  = estInput.menuTab.p_t  ;
p_tk = estInput.menuTab.p_tk ;
t_t  = estInput.menuTab.t_t  ;
t_tk = estInput.menuTab.t_tk ;

% Nodes
ra = estInput.rNodes;
dr = estInput.dNodes;
pb = estInput.bNodes;
nNodes = estInput.nNodes;

% Background consumption
omega = estInput.omega;

% Vector of possible allocations
alphaList = estInput.alphaList';

% Aux
nA = estInput.nA;
nM = estInput.nM;

% Compute indicator function
chosenA = cell(nA,1);
endowment = 100;
for iMenu = 1:nM

    % RDEU of the lottery for each combination of (r,d)

    % Payoff of each allocation
    C1 = ones(nNodes,1)*q_t(iMenu)*endowment*(1-alphaList);
    C2 = ones(nNodes,1)*endowment*alphaList*q_tk(iMenu);

    % Expected utility of each allocation
    EU1 = p_t(iMenu).*( (C1 + omega).^(1 - ra) - omega.^(1 - ra) )./(1 - ra);
    EU2 = p_tk(iMenu).*( (C2 + omega).^(1 - ra) - omega.^(1 - ra) )./(1 - ra);

    % Discounted expected utility of each allocation
    if t_t(iMenu) == 0
        DEU1 = pb .* exp( -dr.*t_t(iMenu) ) .* EU1;
        DEU2 = pb .* exp( -dr.*t_tk(iMenu) ) .* EU2;
    else
        DEU1 = exp( -dr.*t_t(iMenu) ) .* EU1;
        DEU2 = exp( -dr.*t_tk(iMenu) ) .* EU2;
    end
    RDEU  = DEU1 + DEU2;

    % Allocation that maximizes value for each parameter configuration (ra,da)
    [~,idxMax] = max(RDEU,[],2);

    % Indicator matrix defining choice in this menu for each (ra,dr)
    chosenA{iMenu} = sparse((1:nNodes)',idxMax,1,nNodes,nA);

end

% Store what we need for estimation
estInput.chosenA = chosenA;
estInput.nA = nA;

end
