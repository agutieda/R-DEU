function estOutput = estimate_pooled(estInput)
% Estimate parameters of RDEU model, given CTB data, and compute their standard
% errors clustered at the subject level

% Auxiliary objects
% estInput.alphaList = (2.5:5:97.5)'; % Allocation bounds
estInput.alphaList = (0:0.05:1)'; % Allocation midpoints
estInput.nM = size(estInput.menuTab,1);
estInput.nA = length(estInput.alphaList);
estInput.nY = estInput.nA;

% Define integration interval, nodes and weights
rMin = estInput.rMin;
rMax = estInput.rMax;
dMin = estInput.dMin;
dMax = estInput.dMax;
nNodes = estInput.nNodes;

% Halton sequences
ld_sequences = haltonset(2,'Skip',1e3,'Leap',1e2);
ld_sequences = scramble(ld_sequences,'RR2');

% Get the points used for evaluation
ld_sequences = net(ld_sequences,nNodes);

% Distribute integration points uniformly across area
rNodes = ld_sequences(:,1)*(rMax-rMin) + rMin;
dNodes = ld_sequences(:,2)*(dMax-dMin) + dMin;
rNodes(rNodes==1) = 1 + 1e-12;

% Integration weights
intWeights =  ( (rMax-rMin)*(dMax-dMin)/nNodes ).* ones(nNodes,1);

% Store
estInput.rNodes = rNodes;
estInput.dNodes = dNodes;
estInput.intWeights = intWeights;

% Compute Indicator Function
estInput = compute_indicator(estInput);

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
mom.mle.mean_ra    = MOM(1)  ;
mom.mle.median_ra  = MOM(2)  ;
mom.mle.mode_ra    = MOM(3)  ;
mom.mle.sd_ra      = MOM(4)  ;
mom.mle.Z_ra       = MOM(5)  ;
mom.mle.mean_dr    = MOM(6)  ;
mom.mle.median_dr  = MOM(7)  ;
mom.mle.mode_dr    = MOM(8)  ;
mom.mle.sd_dr      = MOM(9)  ;
mom.mle.Z_dr       = MOM(10) ;
mom.mle.mu_ra_plus = MOM(11) ;
mom.mle.p_convex   = MOM(12) ;
mom.mle.iqr_ra     = MOM(13) ;
mom.mle.iqr_dr     = MOM(14) ;


% Use delta method to find std. error. of moments of the estimated distributions
gradient_MOM_hat = jacob_fun(@(x) par2mom(x,estInput),theta_hat);
Cov_MOM_hat= gradient_MOM_hat*Cov_theta_hat*gradient_MOM_hat';

% Std. dev. of estimators
se_MOM_hat   = real( diag( sqrt(Cov_MOM_hat  ) ) )  ;

% Moments
mom.se.mean_ra    = se_MOM_hat(1)  ;
mom.se.median_ra  = se_MOM_hat(2)  ;
mom.se.mode_ra    = se_MOM_hat(3)  ;
mom.se.sd_ra      = se_MOM_hat(4)  ;
mom.se.Z_ra       = se_MOM_hat(5)  ;
mom.se.mean_dr    = se_MOM_hat(6)  ;
mom.se.median_dr  = se_MOM_hat(7)  ;
mom.se.mode_dr    = se_MOM_hat(8)  ;
mom.se.sd_dr      = se_MOM_hat(9)  ;
mom.se.Z_dr       = se_MOM_hat(10) ;
mom.se.mu_ra_plus = se_MOM_hat(11) ;
mom.se.p_convex   = se_MOM_hat(12) ;
mom.se.iqr_ra     = se_MOM_hat(13) ;
mom.se.iqr_dr     = se_MOM_hat(14) ;

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
    round(theta_hat(5),3), ...    
    round(logLike,3), round(pdfCheck,3) );
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