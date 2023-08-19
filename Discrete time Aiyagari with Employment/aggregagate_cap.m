%=========================================================================
%                  Compute Aggregate Capital
%=========================================================================

function [agg_k] = aggregagate_cap(na,agrid,mu)

stationary_k = zeros(na,1);

for ia = 1:na
    id = 2*ia-1;
    stationary_k(ia,1) = mu(id) + mu(id+1);
end

agg_k = stationary_k'*agrid;

end