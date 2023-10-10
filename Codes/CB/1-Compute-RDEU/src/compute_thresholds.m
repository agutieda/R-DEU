function estInput = compute_thresholds(estInput)
% Function to compute the thresholds that make an individual with RDEU
% preferences indifferent between the two alternatives in a risk menu.
% Also computes the discount rate that, conditional on a value of "r", makes the
% individual indifferent between the two alternatives in a time menu.

% Menus
q_t  = estInput.menuTab.q_t  ;
q_tk = estInput.menuTab.q_tk ;
p_t  = estInput.menuTab.p_t  ;
p_tk = estInput.menuTab.p_tk ;
k    = estInput.menuTab.k    ;

% Risk nodes
ra_convex = estInput.ra_convex ;
ra_concave = estInput.ra_concave;

% Background consumption
omega = estInput.omega;

% Vector of possible allocations
alphaList = estInput.alphaList;

% Aux
nA = estInput.nA;
nM = estInput.nM;

% Compute indifference theresholds...
K_convex = nan(estInput.nR_convex,nM);
K_concave  = nan(estInput.nR_concave,nM,nA);
endowment = 100;
for iMenu = 1:nM

    % ... for convex utilities
    C1  = q_t(iMenu)*endowment;
    C2  = q_tk(iMenu)*endowment;
    K_1 = (1/k(iMenu))*log(p_tk(iMenu)/p_t(iMenu));
    num = (omega + C2).^(1-ra_convex) - omega.^(1-ra_convex);
    den = (omega + C1).^(1-ra_convex) - omega.^(1-ra_convex);
    K_2 = (1/k(iMenu)).*log(num./den);
    K_convex(:,iMenu) = K_1 + K_2;

    % ... for concave utilities
    for iA = 1:nA
        alpha = alphaList(iA)/100;
        c_t  = (1-alpha) * endowment * q_t(iMenu);
        c_tk =     alpha * endowment * q_tk(iMenu);
        K_1  =  (1/k(iMenu)) .* log(p_tk(iMenu)/p_t(iMenu)) ;
        K_2  = -(1/k(iMenu)) .* ra_concave.* log( (c_tk+omega)./(c_t+omega) );
        K_3  =  (1/k(iMenu)) .* log(q_tk(iMenu)/q_t(iMenu));
        K_concave(:,iMenu,iA) = K_1 + K_2 + K_3;
    end

end

% Store what we need for estimation
estInput.nM = nM;
estInput.nA = nA;
estInput.K_convex = K_convex;
estInput.K_concave = K_concave;

end
