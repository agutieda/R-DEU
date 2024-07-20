function CovMat_theta_hat = compute_se(loglike, theta_hat, hessian, ...
    cluster_var, se_type, estInput)
% Function to compute the standard errors of NLS estimator

obsTab = estInput.obsTab;

% Number of clusters
N      = size(obsTab,1);
Groups = unique(cluster_var);
M      = length(Groups);
K      = size(hessian,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) Compute J1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get J1 from the Hessian matrix obtained after estimation
J1 = -hessian;

if strcmp(se_type, 'robust')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 2) Compute J2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Set a waitbar
    wait_bar_cluster = waitbar(0,'Computing Robust Std. Errors...',...
        'Name','Computing Std. Errors');
    update_par = 10;
    update_k = round(M/update_par);

    % Gradient of the log-likelihood function for each individual
    J2_i = nan(K,M);
    flag_time=0;
    tic;
    for k = 1:M

        % Individuals in a group
        group_idx = (cluster_var == Groups(k));
        nGroup = sum(group_idx);

        % Data for this individual
        estInput_k = estInput;
        estInput_k.obsTab = obsTab(group_idx,:);

        % Likelihood of indiividuals in that group
        loglike_k = @(x) loglike(x, estInput_k);

        % Jacobian of the log-likelihood evaluated at MLE
        J2_i(:,k) = jacob_fun(loglike_k, theta_hat).*nGroup;

        % Update waitbar
        if mod(k,update_k) == 0
            if flag_time == 0
                time_by_round = toc;
                flag_time = 1;
            end
            update_par = update_par - 1;
            timeLeft = round(time_by_round*update_par/60, 2);
            waitbar(k/M, wait_bar_cluster, ...
                ['Estimated time left: ' num2str(timeLeft),' min']);
        end

    end

    close(wait_bar_cluster);

    % Compute J2 as the sum of the Jacobians of each cluster group
    J2 = zeros(K, K);
    for k=1:M
        J2 = J2 + J2_i(:,k)*J2_i(:,k)';
    end

    J2 = J2./N;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 3) Use J1 and J2 to get the covariance matrix of the se_types
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Asympt_VarMat = J1\J2/J1;

    % Divide by N to get finite sample approximation
    CovMat_theta_hat = Asympt_VarMat./N;

else

    CovMat_theta_hat = inv(J1);

end

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
