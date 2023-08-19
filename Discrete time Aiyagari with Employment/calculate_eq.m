%=========================================================================
%                     Calculate Equilibrium
%=========================================================================

function [r_eq,w_eq,K_eq,dist_eq] = calculate_eq(na,ne,P,agrid,borrow_lim,amax,beta,delta, ...
    alpha,A,L,l_bar)

% tuning parameter
phi = 0.3;

max_demand_excess = 1.0e-8;
min_supply_excess = 1/beta - 1 + delta;

iter = 0; 
cdiff = 1000;
tol_iter = 1.0e-7;
max_iter= 3000;

r_guess = min_supply_excess/2;
r_next = r_guess;

while iter <= max_iter && cdiff>tol_iter
    r = r_next;
    w_r = (1-alpha)*A^(alpha/(1-alpha))*(alpha/(r))^(alpha/(1-alpha));    
    % Aggregate Capital Supply
    [~,cap] = policy_functions(r,w_r,na,borrow_lim,amax);
    mu = stationary_distribution(na,ne,P,agrid,cap);    
    K_s = aggregagate_cap(na,agrid,mu);
    % Aggregate Capital Demand
    K_d = L*l_bar*(alpha*A/r)^(1/(1-alpha));

    % in case there is no convergence, update r
    if K_s > K_d
        min_supply_excess = r;
    else
        max_demand_excess = r;
    end

    cdiff = abs(min_supply_excess-max_demand_excess);

    if mod((iter+1),5) == 0
        disp("Iteration # " + (iter+1) + " - r: " + r + " - Diff: " + cdiff);
    end    

    % update guess
    r_next = phi*min_supply_excess + (1-phi)*max_demand_excess;

    iter = iter +1;

end
r_eq = r;    
w_eq = (1-alpha)*A^(alpha/(1-alpha))*(alpha/(r_eq))^(alpha/(1-alpha));   
K_eq = (K_s + K_d)/2;
dist_eq = mu;

disp("Initial stationary eq.: r = " + r_eq + "; w = " + w_eq + "; K = " + K_eq)
end