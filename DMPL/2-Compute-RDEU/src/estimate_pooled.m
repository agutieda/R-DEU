function estOutput = estimate_pooled(estInput)
% Estimate parameters of RDEU model, given DMPL data, and compute their standard
% errors clustered at the subject level

% Define integration interval, nodes and weights
rMin = -5;
rMax =  5;
estInput.rMin = rMin;
estInput.rMax = rMax;
nNodes = estInput.nNodes;
rNodes = linspace(rMin, rMax, nNodes)';
rNodes(rNodes==1) = 1 + 1e-12;
estInput.rWeights = ones(nNodes, 1) * (rMax - rMin)./ nNodes;
estInput.rNodes = rNodes;
estInput.nNodes = nNodes;

% Compute thresholds for every risk and time menu
estInput = compute_thresholds(estInput);

% Define objective function
obj_fun = @(x) -loglike_fun(x, estInput);

% Optimize over a list of initial values using fminsearch
theta_hat_i_list = nan(size(estInput.theta_0));
loglike_i_list   = nan(size(estInput.theta_0, 1), 1);
for i = 1:size(estInput.theta_0, 1)

    [theta_hat_i_list(i, :), loglike_i_list(i)] = ...
        fminsearch(obj_fun, estInput.theta_0(i, :), estInput.opt_x0);

    [~, ~, ~, pdfCheck_i] = loglike_fun(theta_hat_i_list(i, :), estInput);

    if abs(pdfCheck_i-1)>0.01
        loglike_i_list(i) = nan;
    end

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

% Loglikelihood at MLE
[logLike, logLike_R, logLike_T, pdfCheck, rho_L1_R, rho_L1_T, ...
    rhoObs_L1_R, rhoObs_L1_T] = loglike_fun(theta_hat, estInput);

% If discretized PDF doesn't add to 1, refine estimation
if abs(pdfCheck-1) > 0.01 || theta_hat(2) < 0.05

    % Increase number of nodes and center around current estimates
    rMin   = theta_hat(1) - 10*theta_hat(2);
    rMax   = theta_hat(1) + 10*theta_hat(2);
    nNodes = estInput.nNodes;
    rNodes = linspace(rMin, rMax, nNodes)';
    rNodes(rNodes==1) = 1 + 1e-12;
    estInput.rWeights = ones(nNodes, 1) * (rMax - rMin)./ nNodes;
    estInput.rNodes = rNodes;
    estInput.nNodes = nNodes;

    % Compute thresholds for every risk and time menu again
    estInput = compute_thresholds(estInput);

    % Define objective function again
    obj_fun = @(x) -loglike_fun(x, estInput);

    % Maximize log-likelihood function using a quasi-newton method
    [theta_hat, ~, exitflag, ~, ~, hessian_loglike] = ...
        fminunc(obj_fun, theta_hat, estInput.optFminunc_mle);

    [logLike, logLike_R, logLike_T, pdfCheck, rho_L1_R, rho_L1_T, ...
        rhoObs_L1_R, rhoObs_L1_T, ] = loglike_fun(theta_hat, estInput);

    if abs(pdfCheck-1) > 0.01
        warning( ...
            ['Integration is not adding to 1 in pooled estimates: ',...
            'Revise bounds and number of nodes of the support of r .'] );
    end

end

% Compute robust standard errors of estimated parameters
cluster_var   = estInput.obsTab.subjectID;
Cov_theta_hat = compute_se(@loglike_fun, ...
    theta_hat, -hessian_loglike, cluster_var, 'robust', estInput);
se_theta_hat = diag(sqrt(Cov_theta_hat))' ;

if size(theta_hat,2) == 4
    theta_hat = [theta_hat,0];
    se_theta_hat = [se_theta_hat,nan];
end

% Display results
fprintf('Pooled Estimates:\n');
fprintf(...
    'mu_ra:%g sig_ra:%g mu_dr:%g sig_dr:%g rho:%g logLike:%g pdfCheck:%g \n', ...
    round(theta_hat(1),3), round(theta_hat(2),3), ...
    round(theta_hat(3),3), round(theta_hat(4),3), ...
    round(theta_hat(5),3), round(logLike,3), round(pdfCheck,3) );
fprintf('\n');

% Store results in structure
estOutput = struct( ...
    'theta_hat', theta_hat, ...
    'se_theta_hat', se_theta_hat, ...
    'logLike', logLike, ...
    'logLike_R', logLike_R, ...
    'logLike_T', logLike_T, ...
    'rho_L1_R', rho_L1_R, ...
    'rho_L1_T', rho_L1_T, ...
    'rhoObs_L1_R', rhoObs_L1_R, ...
    'rhoObs_L1_T', rhoObs_L1_T, ...
    'exitFlag', exitflag, ...
    'nObs', size(estInput.obsTab,1), ...
    'pdfCheck', pdfCheck);

end