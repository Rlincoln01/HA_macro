%=========================================================================
%                  Compute Stationary distribution
%=========================================================================

function [mu,T] = stationary_distribution(na,ne,P,agrid,capital)

% compute the size of the grid
ns = na*ne; 

% Matrix that indicates, given the current state, capital chosen
% in the next period
next_a_index = zeros(na,ne); 
for ia = 1:na
    for ie=1:ne
        s = capital(ia,ie);
        [~,id] = min((agrid-s).^2);
        next_a_index(ia,ie) = id;
    end
end
% T is the transition matrix, which is a a markov chain transition matrix
T = zeros(ns,ns);


for ia=1:na
    for ie = 1:ne
        for ie_next = 1:ne
            T(2*(ia-1) + ie, 2*(next_a_index(ia,ie)-1) + ie_next) = P(ie,ie_next);
        end
    end
end   

mu = asymptotics(dtmc(T));
end