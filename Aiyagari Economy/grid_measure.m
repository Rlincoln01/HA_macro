% ========================================================================
%                   Income Process Generator Matrix
% ========================================================================
% 
% Calculates the measure for grid, used for integration purposes  
%
% 
% ========================================================================

function [dm_tilde,dmf,dmb] = grid_measure(grid,n)
    dmf = ones(n,1);                      % measure forward difference
    dmb = ones(n,1);                      % measure backward difference
    dmf(1:n-1) = grid(2:n) - grid(1:n-1); % vector of grid-varying measures (F)
    dmb(2:n) = grid(2:n) - grid(1:n-1);   % vector of grid-varying measures (B)
    dmf(n) = dmf(n-1); dmb(1) = dmb(2);

    dm_tilde = 0.5.*(dmf+dmb);
    dm_tilde(1) = 0.5*dmf(1); dm_tilde(n) = 0.5*dmb(n);
end