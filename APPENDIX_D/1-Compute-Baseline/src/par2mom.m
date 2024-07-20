function Y = par2mom(X,estInput)
% Function that computes the mean, median, mode and standard deviation of a
% truncated normal distribution, given values of scale and location parameters

mu_ra  = X(1);
sig_ra = X(2);
mu_dr  = X(3);
sig_dr = X(4);

% Mean of distribution of risk aversion, conditional on ra>0
alf = (0-mu_ra)/sig_ra;
mu_ra_plus = mu_ra + sig_ra*normpdf(alf)/(1 - normcdf(alf));

% IQR: ra
aux_fun = @(x) normcdf(x,mu_ra,sig_ra);
p25 = fzero(@(x) aux_fun(x) - 0.25, [estInput.rMin,estInput.rMax]);
p75 = fzero(@(x) aux_fun(x) - 0.75, [estInput.rMin,estInput.rMax]);
iqr_ra = p75 - p25;

% IQR: dr
aux_fun = @(x) normcdf(x,mu_dr,sig_dr);
p25 = fzero(@(x) aux_fun(x) - 0.25, [-10,10]);
p75 = fzero(@(x) aux_fun(x) - 0.75, [-10,10]);
iqr_dr = p75 - p25;

% Output
Y = [mu_ra_plus; iqr_ra; iqr_dr];

end
