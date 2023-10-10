function [mean_x,median_x,mode_x,std_x,Z_x] = mom_dr(mu,sigma,a,b)
    % Function to compute some of the moments of ra
    %
    % This file: Truncated Normal Distribution
    % The distribution is characterized by a location parameter "mu", a scale
    % parameter "sigma", lower bound "a" and upper bound "b"
    
    % Uses a procedure to compute the mean and std. dev. suggested by Jorge Fernandez-de-Cossio-Diaz
    % which avoids catastrophic cancellation when mu lies far away from the integration bounds
    
    % Auxiliary parameters
    alpha = (a-mu)./sigma ;
    beta  = (b-mu)./sigma ;
    Z_x = normcdf(beta) - normcdf(alpha);
    
    % Mean
    mean_x = mu + sigma*sqrt(2/pi)*F_1(alpha/sqrt(2),beta/sqrt(2));
    
    % Std. Deviation
    std_x = sigma * sqrt( 1 + (2/sqrt(pi))*F_2(alpha/sqrt(2),beta/sqrt(2)) ...
        - (2/pi)*( F_1(alpha/sqrt(2),beta/sqrt(2))^2 ) );
    
    % Median
    aux = min((normcdf(alpha)+normcdf(beta))./2 , 1-eps);
    median_x = mu + sigma.*norminv( aux );
    
    % Mode
    mode_x = a.*(mu<a) + b.*(mu>b) + mu.*(a<=mu && mu<=b);
    
    end
    
    %% Auxiliary functions
    
    function f = F_1(x,y,varargin)
    
    if ~isempty(varargin)
        thereshold = varargin{1};
    else
        thereshold = 1e-7;
    end
    
    if abs(x) > abs(y)
        aux = x;
        x = y;
        y = aux;
    elseif isinf(y)
        f = sign(y)/ erfcx(sign(y)*x);
        return;
    elseif abs(x-y) <= thereshold
        epsilon = y-x;
        f = sqrt(pi)*x + (sqrt(pi)/2 + (-sqrt(pi)*x/6 + (-sqrt(pi)/12 + x*(sqrt(pi)/90 + ...
            (sqrt(pi)*x^2)/90)*epsilon)*epsilon)*epsilon)*epsilon ;
        return;
    end
    
    Delta = exp(x^2-y^2);
    
    if max(x,y) < 0
        f = (1-Delta) / (Delta*erfcx(-y) - erfcx(-x));
    elseif min(x,y)>0 || y==Inf
        f = (1-Delta) / ( erfcx(x) - Delta*erfcx(y) );
    else
        f = exp(-x^2) * (1-Delta) / (erf(y)-erf(x));
    end
    
    end
    
    
    function f = F_2(x,y,varargin)
    
    if ~isempty(varargin)
        thereshold = varargin{1};
    else
        thereshold = 1e-7;
    end
    
    if abs(x) > abs(y)
        aux = x;
        x = y;
        y = aux;
    elseif (x==Inf && y==-Inf) || (x==-Inf && y==Inf)
        f = 0;
        return;
    elseif isinf(y)
        f = sign(y)*x / erfcx(sign(y)*x);
        return;
    elseif abs(x-y) <= thereshold
        epsilon = y-x;
        f = sqrt(pi)*x^2 - sqrt(pi)/2 + (sqrt(pi)*x + (sqrt(pi)/3 - sqrt(pi)*x^2/3 +  ...
            (((sqrt(pi)/30 + sqrt(pi)*x^2/45)*x^2 - 4*sqrt(pi)/45)*epsilon - ...
            sqrt(pi)*x/3)*epsilon)*epsilon)*epsilon;
        return;
    end
    
    Delta = exp(x^2-y^2);
    
    if max(x,y) < 0
        f = (x - Delta*y) / (Delta*erfcx(-y) - erfcx(-x));
    elseif min(x,y)>0
        f = (x - Delta*y) / ( erfcx(x) - Delta*erfcx(y) );
    else
        f = exp(-x^2) * ( x-Delta*y ) / (erf(y)-erf(x));
    end
    
    end
    
    