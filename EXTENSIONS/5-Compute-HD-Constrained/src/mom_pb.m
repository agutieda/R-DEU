function [mean_x,median_x,mode_x,std_x,Z_x] = mom_pb(alf,bet,a,b)
    % Function to compute some of the moments of pb
    %
    % This file: Truncated Beta Distribution
    % The distribution is characterized by shape parameters "alf" and "bet",
    % lower bound "a" and upper bound "b"
    
    % Make sure parameters are positive
    if alf<=0 
        alf = 1e-8;
    end
    
    if bet<=0
        bet = 1e-8;
    end
    
    Z_x = betacdf(b,alf,bet)-betacdf(a,alf,bet);
    
    % Mean and Std. Deviation
    [mean_x,var_x] = betastat(alf,bet);
    std_x = sqrt(var_x);
    
    % Median
    if alf>1 && bet>1
        median_x = ( alf-(1/3) )/( alf+bet-(2/3) ) ;
    elseif alf==bet
        median_x = 0.5;
    elseif alf<1 && bet<1
        median_x = nan;
    elseif alf<0.01 || bet > 1e4
        median_x = 0;
    elseif  bet<0.01 || alf > 1e4
        median_x = 1;
    else    
        % Brute force
        auxFun = @(x) betainc(x,alf,bet)-0.5;
        options = optimset('Display','off'); 
        median_x = fzero(auxFun,[0,1],options);   
    end
    
    % Mode
    if alf>1 && bet>1
        mode_x = (alf-1) / (alf+bet-2);
    elseif alf==1 && bet==1
        mode_x = 0.5; % Any value in (0,1)
    elseif alf<=1 && bet>1
        mode_x = 0;
    elseif alf>1 && bet<=1
        mode_x = 1;
    else
        mode_x = nan;
    end
    
end
