function estOutput = estimate_pooled(estInput)
% Estimate parameters of model using NLS, given CTB data, and compute their standard
% errors clustered at the subject level

% Auxiliary objects
estInput.nM = size(estInput.menuTab,1);

% Define objective function
obj_fun = @(x) sse_fun(x, estInput);

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

% Predicted choice
[~,a_hat] = sse_fun(theta_hat, estInput);

% Warn in case of non-convergence
if exitflag <= 0
    warning('MLE did not converge');
end

% Compute robust standard errors of estimated parameters
cluster_var   = estInput.obsTab.subjectID;
Cov_theta_hat = compute_se(@sse_fun, ...
    theta_hat, -hessian_loglike, cluster_var, 'robust', estInput);
se_theta_hat = diag(sqrt(Cov_theta_hat))' ;

% Display results
fprintf('Pooled Estimates:\n');
fprintf(...
    'mu_ra:%g mu_dr:%g \n', ...
    round(theta_hat(1),3), round(theta_hat(2),3));
fprintf('\n');

% Store results in structure
estOutput = struct( ...
    'theta_hat', theta_hat, ...
    'se_theta_hat', se_theta_hat, ...
    'exitFlag', exitflag, ...
    'a_hat', a_hat, ...
    'nObs', size(estInput.obsTab,1) );

end
