function Y = par2mom(X,estInput)
% Function that computes the mean, median, mode and standard deviation of a
% truncated normal distribution, given values of scale and location parameters

rMax = estInput.rMax;
rMin = estInput.rMin;
dMax = estInput.dMax;
dMin = estInput.dMin;
bMax = estInput.bMax;
bMin = estInput.bMin;


mu_ra  = X(1);
sig_ra = X(2);
mu_dr  = X(3);
sig_dr = X(4);
mu_pb  = X(5);
sig_pb = X(6);

% Moments of distribution of risk aversion
[mean_ra,median_ra,mode_ra,sd_ra,Z_ra] = mom_ra(mu_ra,sig_ra,rMin,rMax);

% Moments of distribution of discount rate
[mean_dr,median_dr,mode_dr,sd_dr,Z_dr] = mom_dr(mu_dr,sig_dr,dMin,dMax);

% Moments of distribution of present bias
[mean_pb,median_pb,mode_pb,sd_pb,Z_pb] = mom_pb(mu_pb,sig_pb,bMin,bMax);

% Mean of distribution of risk aversion, conditional on ra>0
alf = (0-mu_ra)/sig_ra;
mu_ra_plus = mu_ra + sig_ra*normpdf(alf)/(1 - normcdf(alf));

% Probability of convex utility
p_convex = 100*normcdf(0, mu_ra, sig_ra);

% IQR: ra
aux_fun = @(x) cdf_ra(x,mu_ra,sig_ra,rMin,rMax);
p25 = fzero(@(x) aux_fun(x) - 0.25, [rMin,rMax]);
p75 = fzero(@(x) aux_fun(x) - 0.75, [rMin,rMax]);
iqr_ra = p75 - p25;

% IQR: dr
aux_fun = @(x) cdf_dr(x,mu_dr,sig_dr,dMin,dMax);
p25 = fzero(@(x) aux_fun(x) - 0.25, [dMin,dMax]);
p75 = fzero(@(x) aux_fun(x) - 0.75, [dMin,dMax]);
iqr_dr = p75 - p25;

% IQR: pb
aux_fun = @(x) cdf_pb(x,mu_pb,sig_pb,bMin,bMax);
p25 = fzero(@(x) aux_fun(x) - 0.25, [bMin,bMax]);
p75 = fzero(@(x) aux_fun(x) - 0.75, [bMin,bMax]);
iqr_pb = p75 - p25;


% Output
Y = [
    mean_ra; median_ra; mode_ra; sd_ra; Z_ra;
    mean_dr; median_dr; mode_dr; sd_dr; Z_dr;
    mean_pb; median_pb; mode_pb; sd_pb; Z_pb;
    mu_ra_plus; p_convex; iqr_ra; iqr_dr; iqr_pb;
    ];


end
