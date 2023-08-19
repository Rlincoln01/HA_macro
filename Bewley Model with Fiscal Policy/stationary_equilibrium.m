%=========================================================================
%                   Solve for the Stationary Equilibrium
%=========================================================================

%=========================================================================
% Solves the system: B(r) = S(r)
% Returns:
% - r
% - B
%=========================================================================


function [B_eq,r_eq] = stationary_equilibrium(ydist)

%=========================================================================
%                   Parameters and technicalities
%=========================================================================

% == Household Parameters == %
rho = .05; % beta ~ 0.99
y_1 = 0.1; y_2 = 0.15; y_3 = 0.2;
y = [y_1,y_2,y_3];
tau_0 = 0.2;

% == Iteration parameters == %
iter = 0; 
cdiff = 1000;
tol_iter = 1.0e-7;
max_iter= 3000;
phi = 0.9; % tuning parameter

% == r guess & update == %
min_supply_excess = 1.0e-8;
max_demand_excess = rho;
rguess = max_demand_excess/2;
r_next = rguess; % start loop

% == Asset Grid == %
I=500; % size of asset grid
borrow_lim = -0.15;
amin = borrow_lim;
amax = 3;
a = linspace(amin,amax,I)'; % asset grid
da = (amax-amin)/(I-1); % step on asset grid

%=========================================================================
%                   Equilibrium Loop
%=========================================================================


while iter <= max_iter && cdiff>tol_iter
    r= r_next;

    % == Aggregate bond demand == %
    [~,~,~,A] = HJB_solve_new(r);
    [~,g] = KF_solve_new(A);
    B_d = g(:,1)'*a*da + g(:,2)'*a*da + g(:,3)'*a*da;

    % == Aggregate Bond Supply == %
    taxes = tau_0*(y*ydist);
    B_s = taxes/r;
    
    % == If no convergence, update r == %
    if B_d > B_s
       max_demand_excess = r;
    else
       min_supply_excess = r;
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
B_eq = (B_s + B_d)/2;

disp("Initial stationary eq.: r = " + r_eq + " ; B = " + B_eq)

end





