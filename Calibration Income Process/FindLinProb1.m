%=========================================================================
% FindLinProb1
%   This function takes a value outine takes in a value xi and an array x,
% and finds the two indices in x that correspond to the values immediately
% to the left and right of xi. It then returns the two indices as y, and
% the associated probabilities of the value being between these two indices
% as p.
%
%
% Author: Greg Kaplan, G. Violante and Ben Moll
%=========================================================================
function [y, p] = FindLinProb1(xi, x)

n = length(x);
if n == 1
    p = [1.0, 0.0];
    y = [1, 1];
    return;
end

[~, locL] = max(xi > x);
if xi <= x(1)
    y = [1, 2];
    p = [1.0, 0.0];
elseif locL >= n
    locL = n - 1;
    y = [n - 1, n];
    p = [0.0, 1.0];
elseif x(locL + 1) == x(locL)
    y = [locL, locL + 1];
    p = [0.5, 0.5];
else
    y = [locL, locL + 1];
    p(2) = (xi - x(locL)) / (x(locL + 1) - x(locL));
    p(1) = 1.0 - p(2);
end

end
