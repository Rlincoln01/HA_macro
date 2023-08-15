%=========================================================================
% SymmetricPowerSpacedGrid
%
% Gives a grid spaced between center-width and center+width based on the 
% interval [-1,1] with a function x^(1/k) either side. If k = 1, the grid
% is linear; and if k = 0, is L-shaped
%
% Obs: n>= 2
% Author: Greg Kaplan, G. Violante and Ben Moll
%=========================================================================

function[y] = SymmetricPowerSpacedGrid(n,k,width,center)

if n == 2
    y(1) = center - width;
    y(2) = center + width;
end

% create a linear spaced grid between -1,1
x = linspace(-1,1,n);

for i=1:n
    if x(i)>0
        z(i) = x(i)^(1/k);
    elseif x(i) == 0
        z(i) = 0;
    else 
        z(i) = -((-x(i))^(1/k));
    end
y= center + width.*z;    
end