function grid = power_grid(par_min,par_max,n)
    x = linspace(0,1,n)';           % equispaced grid on [0,1]
    coeff = 50;                     % controls the density of grid in lower part
    power = 5;                      % controls the density of grid in upper part
    
    % x grid
    xx = x + coeff*x.^power;        
    xmax = max(xx); xmin = min(xx); 
    
    grid = (par_max-par_min)/(xmax-xmin)*xx + par_min; 
end
