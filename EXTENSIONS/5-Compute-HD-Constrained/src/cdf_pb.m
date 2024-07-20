function cdf_x = cdf_pb(x,mu,sigma,a,b)
% Function to evaluate the cdf of pb
%
% This file: Truncated Beta Distribution
% The distribution is characterized by shape parameters "alpha" and "beta",
% lower bound "a" and upper bound "b"

xi    = (x-mu)./sigma;
alpha = (a-mu)./sigma;
beta  = (b-mu)./sigma;
Z = normcdf(beta)-normcdf(alpha);
cdf_x = (normcdf(xi)-normcdf(alpha))./Z;
cdf_x(x<a)=0;
cdf_x(x>b)=1;

% Z = betacdf(b,alpha,beta)-betacdf(a,alpha,beta);
% cdf_x = (betacdf(x,alpha,beta)-betacdf(a,alpha,beta))./Z;
% cdf_x(x<a)=0;
% cdf_x(x>b)=1;

end