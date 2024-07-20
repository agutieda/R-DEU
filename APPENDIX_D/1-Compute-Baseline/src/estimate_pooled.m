function estOutput = estimate_pooled(estInput)
% Estimate parameters of RDEU model, given CTB data, and compute their standard
% errors clustered at the subject level

% Auxiliary objects
estInput.alphaList = (2.5:5:97.5)'; % Allocation bounds
estInput.nM = size(estInput.menuTab,1);
estInput.nA = length(estInput.alphaList);
estInput.nY = estInput.nA + 1;

% Define integration interval, nodes and weights
rMin = estInput.rMin;
rMax = estInput.rMax;
nNodes = estInput.nNodes;
rNodes = linspace(rMin, rMax, nNodes)';
rNodes(rNodes==1) = 1 + 1e-12;

% Separate convex and concave part
ra_convex = rNodes(rNodes<=0);
nR_convex = length(ra_convex);
ra_concave = rNodes(rNodes>0);
nR_concave = length(ra_concave);

% Integration weights
intWeights_convex  =  ( (0-rMin)/nR_convex ).* ones(nR_convex,1);
intWeights_concave =  ( (rMax-0)/nR_concave ).* ones(nR_concave,1);

% Store
estInput.ra_convex = ra_convex;
estInput.ra_concave = ra_concave;
estInput.intWeights_convex = intWeights_convex;
estInput.intWeights_concave = intWeights_concave;
estInput.nR_convex = nR_convex;
estInput.nR_concave = nR_concave;

% Compute thresholds for every menu
estInput = compute_thresholds(estInput);

% Define objective function
obj_fun = @(x) -loglike_fun(x, estInput);

% Optimize over a list of initial values using fminsearch
theta_hat_i_list = nan(size(estInput.theta_0));
loglike_i_list   = nan(size(estInput.theta_0, 1), 1);
for i = 1:size(estInput.theta_0, 1)
    [theta_hat_i_list(i, :), loglike_i_list(i)] = ...
        fminsearch(obj_fun, estInput.theta_0(i, :), estInput.opt_x0);
end

% Choose value with lowest log-likelihood
[~, i_min] = min(loglike_i_list);
theta_0 = theta_hat_i_list(i_min, :);

% Maximize log-likelihood function using a quasi-newton method
[theta_hat, ~, exitflag, ~, ~, hessian_loglike] = ...
    fminunc(obj_fun, theta_0, estInput.optFminunc_mle);

% Warn in case of non-convergence
if exitflag <= 0
    warning('MLE did not converge');
end

% Compute robust standard errors of estimated parameters
cluster_var   = estInput.obsTab.subjectID;
Cov_theta_hat = compute_se(@loglike_fun, ...
    theta_hat, -hessian_loglike, cluster_var, 'robust', estInput);
se_theta_hat = diag(sqrt(Cov_theta_hat))' ;

% Moments of the estimated distributions
MOM = par2mom(theta_hat,estInput);
mom.mle.mu_ra_plus = MOM(1) ;
mom.mle.iqr_ra     = MOM(2) ;
mom.mle.iqr_dr     = MOM(3) ;

% Use delta method to find std. error. of moments of the estimated distributions
gradient_MOM_hat = jacob_fun(@(x) par2mom(x,estInput),theta_hat);
Cov_MOM_hat= gradient_MOM_hat*Cov_theta_hat*gradient_MOM_hat';

% Std. dev. of estimators
se_MOM_hat   = real( diag( sqrt(Cov_MOM_hat  ) ) )  ;

% Moments
mom.se.mu_ra_plus = se_MOM_hat(1) ;
mom.se.iqr_ra     = se_MOM_hat(2) ;
mom.se.iqr_dr     = se_MOM_hat(3) ;

% Loglikelihood at MLE
[logLike, rhoY_hat, rhoY_obs, pdfCheck] = loglike_fun(theta_hat, estInput);

if abs(pdfCheck-1) > 0.01
    warning( ...
        ['Integration is not adding to 1 in pooled estimates: ',...
        'Revise bounds and number of nodes of the support of r .'] );
end

% Display results
fprintf('Pooled Estimates:\n');
fprintf(...
    'mu_ra:%g sig_ra:%g mu_da:%g sig_da:%g rho:%g logLike:%g pdfCheck:%g \n', ...
    round(theta_hat(1),3), round(theta_hat(2),3), ...
    round(theta_hat(3),3), round(theta_hat(4),3), ...
    round(theta_hat(5),3), round(logLike,3), round(pdfCheck,3) );
fprintf('\n');

% Store results in structure
estOutput = struct( ...
    'theta_hat', theta_hat, ...
    'se_theta_hat', se_theta_hat, ...
    'mom', mom, ...
    'logLike', logLike, ...
    'rhoY_hat', rhoY_hat, ...
    'rhoY_obs', rhoY_obs, ...
    'exitFlag', exitflag, ...
    'nObs', size(estInput.obsTab,1), ...
    'pdfCheck', pdfCheck);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Auxiliary function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function jac = jacob_fun(f, x)
% Compute the Jacobian numerically using forward differences
step = 1e-5;
step_inv = 1/step;
nx = length(x); % Dimension of the input x;
f0 = feval(f, x); % caclulate f0, when no perturbation happens
jac = zeros(length(f0), nx);
% Do perturbation
for i = 1 : nx
    xplus = x;
    xplus(i) =  x(i) + step;
    jac(:, i) = (feval(f, xplus) - f0) .* step_inv;
end
end