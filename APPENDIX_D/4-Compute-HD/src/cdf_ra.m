function cdf_x = cdf_ra(x,mu,sigma,a,b)
% Function to evaluate the cdf of ra
%
% This file: Truncated Normal Distribution
% The distribution is characterized by a location parameter "mu", a scale
% parameter "sigma", lower bound "a" and upper bound "b"

xi    = (x-mu)./sigma;
alpha = (a-mu)./sigma;
beta  = (b-mu)./sigma;
Z = normcdf(beta)-normcdf(alpha);
cdf_x = (normcdf(xi)-normcdf(alpha))./Z;
cdf_x(x<a)=0;
cdf_x(x>b)=1;

end