%=========================================================================
% Compute ergodic distribution of a CTMC
%
% Inputs: 
% 1)TM -  Transition Matrix
% 2)n  -  Number of grid points
%
% Output: g_tilde ~ Stationary Distribution
%
% Author: Rafael Lincoln 
%=========================================================================

function [g_tilde] = stationary_dist(TM,n)

% take the transpose of the infinitesimal generator
AT = TM';

% fix one value of the matrix, otherwise it has no inverse (its singular)
% i_fix = (n+1)/2;
i_fix = 1;
row = [zeros(1,i_fix-1),1,zeros(1,n-i_fix)];
AT(i_fix,:) = row;

% Set the array of zeroes and also fix a value
b = zeros(n,1);
b(i_fix) = 1;

% solve for the stationary distribution
g_tilde = AT\b;
g_sum = g_tilde'*ones(n,1);
g_tilde = g_tilde./g_sum;


end