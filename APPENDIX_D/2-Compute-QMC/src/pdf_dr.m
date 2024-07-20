function pdf_x = pdf_dr(x,mu,sigma,a,b)
    % Function to evaluate the pdf of dr
    %
    % This file: Truncated Normal Distribution
    % The distribution is characterized by a location parameter "mu", a scale
    % parameter "sigma", lower bound "a" and upper bound "b"
    
    xi    = (x-mu)./sigma ;
    alpha = (a-mu)./sigma ;
    beta  = (b-mu)./sigma ;
    Z = normcdf(beta)-normcdf(alpha);
    pdf_x = normpdf(xi) ./ (sigma*Z);
    pdf_x(x<a | x>b)=0;
    
    
    end