function pdf_x = pdf_pb(x,mu,sigma,a,b)
    % Function to evaluate the pdf of pb
    %
    % This file: Truncated Beta Distribution
    % The distribution is characterized by shape parameters "alpha" and "beta",
    % lower bound "a" and upper bound "b"
    
    xi    = (x-mu)./sigma ;
    alpha = (a-mu)./sigma ;
    beta  = (b-mu)./sigma ;
    Z = normcdf(beta)-normcdf(alpha);
    pdf_x = normpdf(xi) ./ (sigma*Z);
    pdf_x(x<a | x>b)=0;
    
    % Z = betacdf(b,alpha,beta)-betacdf(a,alpha,beta);
    % pdf_x = betapdf(x,alpha,beta) ./ Z;
    % pdf_x(x<a | x>b)=0;
    
    
end