%=;========================================================================
%                     Transition Prices
%=========================================================================

function [r, w] = transition_prices(K,A,alpha,L,l_bar)

r = A.*alpha.*(K./l_bar.*L).^(alpha-1);

w = A.*(1-alpha).*(K./l_bar.*L).^(alpha);

end