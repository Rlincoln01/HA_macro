%=========================================================================
%                   Stationary Equilibrium
%=========================================================================
% Description: Solves the recursive equilibrium of the Aiyagari economy
% 
%
% Inputs: 
%
% - Calibration of the income process:
%   • Permanent Income parameters
%   • Transitory Income parameters
% - Asset grid:
%   • Object which defines the support of the distribution g
% - Model Specification:
%   • Open Economy
%   • Individual Fixed affects
%   • Discount rate heterogeneity 
%
% Outputs: Equilibrium objects of the economy
% - Aggregate variables: Y,C,S,K,L,lump_sum
% - Prices: r,w
% - Densities: g, a (marginal wealth)
%
%=========================================================================

function ss = stationary_equilibrium(calibration,parameters,specification,a_grid)

% options
display_iterations = 0;

%% ============== 0 step) Define Parameters =========================== %%

rho = calibration(1);
delta_rho = calibration(2);
tau = calibration(3);
borrowing_limit = calibration(4);

% 1.1 - Calibration/Parametrization
alpha = parameters.alpha;
r_open_econ = parameters.r_open_econ;
dep = parameters.dep;
n_p = parameters.n_p;
n_t = parameters.n_t;
n_fe = parameters.n_fe;
sigma_omega = parameters.sigma_omega;
n_rho = parameters.n_rho;


% Asset grid
n_a = length(a_grid);
da_tilde = grid_measure(a_grid,n_a); % measure

%% ============== 1st) Labor Market Equilibrium ======================= %%

% 1.1 - Obtain ergodic distribution of productivity
n_z = n_p*n_t;
[lambda_z,zgrid] = inc_matrix(parameters);
dz_tilde = grid_measure(zgrid,n_z); % Define measure for the productivity grid for integration
z_dist = stationary_dist(lambda_z,n_z,dz_tilde);

% 1.2 - obtain discrete distribution of FE, if applies
if specification.ind_fixed_effects == 1
    % Discrete dist + grid + measure
    [omega_dist,omega_grid,domega_tilde] = FE_grid(n_fe,sigma_omega); 
end

% 1.3 - Labor Market clearing
if specification.ind_fixed_effects == 0
    Labor = sum(exp(zgrid).*z_dist.*dz_tilde);
else
    Labor = sum(exp(omega_grid').*omega_dist.*domega_tilde)*sum(exp(zgrid).*z_dist.*dz_tilde);
end

% Do bisection method if closed economy
%=========================================================================
%                   Equilibrium Loop
%=========================================================================

if specification.closed_economy == 1
    % == Iteration parameters == %
    iter = 0; 
    cdiff = 1000;
    tol_iter = 1.0e-7;
    max_iter= 3000;
    phi = 0.5; % tuning parameter
    min_supply_excess = -dep;
    max_demand_excess = rho;
    rguess = (max_demand_excess + min_supply_excess)/2;
    r_next = rguess; % initial guess is the open economy rate
    while iter <= max_iter && cdiff >tol_iter
        r_iter=r_next;

        % =========== Calculate Equilibrium wage =========== %
        w = (1-alpha)*(alpha/(r_iter+dep))^(alpha/(1-alpha));
        
        % =========== Calculate lump_sum transfers =========== %
        lump_sum = tau*w*Labor;

        % =========== Household Block =========== %
       
        if specification.disc_rate_heterogeneity == 1
            % preallocation for policy functions and distributions
            g_dist = zeros(n_a,n_z,n_rho);
            con_dist = zeros(n_a,n_z,n_rho);
            for j=1:n_rho
                rho_j = rho + (j-3)*delta_rho;
                if specification.ind_fixed_effects == 0
                    [g_dist(:,:,j),con_dist(:,:,j)] = HJB(parameters,0,a_grid, ...
                        borrowing_limit,r_iter, ...
                        w,lump_sum,tau,rho_j);
                else
                    g_cell = cell(n_fe);
                    con_cell = cell(n_fe);
                    sav_cell = cell(n_fe);
                    % Loop for fixed effects (solve for each)
                    for i=1:n_fe
                        [x,y,z] = HJB(parameters,omega_grid(i), ...
                            a_grid,borrowing_limit,r_iter,w,lump_sum,tau,rho_j);
                        g_cell{i} = x;
                        con_cell{i} = y;
                        sav_cell{i} = z;
                    end
                    g_final = 0;
                    con_final = 0;
                    sav_final = 0;
                    % calculate aggregate consumption, capital and distribution
                    for i = 1:n_fe
                        g_final = g_final + omega_dist(i).*g_cell{i}.*domega_tilde(i);
                        con_final = con_final + omega_dist(i).*con_cell{i}.*domega_tilde(i);
                        sav_final = sav_final + omega_dist(i).*sav_cell{i}.*domega_tilde(i);
                    end
                    g_dist(:,:,j) = g_final;
                    con_dist(:,:,j) = con_final;
                end
            end
            con_final = mean(con_dist,3);
            g_final =mean(g_dist,3);
        else
            if specification.ind_fixed_effects == 0
                [g_final,con_final,~] = HJB(parameters,0,a_grid, ...
                    borrowing_limit,r_iter,w,lump_sum,tau);
            else
                g_cell = cell(n_fe);
                con_cell = cell(n_fe);
                sav_cell = cell(n_fe);
                % Loop for fixed effects (solve for each)
                for i=1:n_fe
                    [x,y,z] = HJB(parameters,omega_grid(i), ...
                        a_grid,borrowing_limit,r_iter,w,lump_sum,tau);
                    g_cell{i} = x;
                    con_cell{i} = y;
                    sav_cell{i} = z;
                end
                g_final = 0;
                con_final = 0;
                sav_final = 0;
                % calculate aggregate consumption, capital and distribution
                for i = 1:n_fe
                    g_final = g_final + omega_dist(i).*g_cell{i}.*domega_tilde(i);
                    con_final = con_final + omega_dist(i).*con_cell{i}.*domega_tilde(i);
                    sav_final = sav_final + omega_dist(i).*sav_cell{i}.*domega_tilde(i);
                end
            end
        end    
        % =========== Calculate Aggregate demand and supply of capital =========== %
        
        % Capital Supply
        a_marg_dist = g_final*dz_tilde; % marginal wealth dist
        K_s = sum(a_grid.*a_marg_dist.*da_tilde);

        % Capital Demand
        K_d = Labor*(alpha/(r_iter+dep))^(1/(1-alpha)); % check to see if

        % == If no convergence,update == %
        if K_s > K_d
            max_demand_excess = r_iter;
        else
            min_supply_excess = r_iter;
        end

        cdiff = abs(min_supply_excess-max_demand_excess);
        K_diff = K_d - K_s;
        if (mod((iter+1),5) == 0) && display_iterations == 1
            disp("Iteration # " + (iter+1) + " - r: " + r_iter + " - Diff: " + cdiff + " given Kdiff: " + K_diff);
        end    
    
        % update guess
        r_next = phi*min_supply_excess + (1-phi)*max_demand_excess;
    
        iter = iter +1;
    end
    % Equilibrium interest rate
    r_eq = r_iter;
    % Equilibrium assets
    K_eq = (K_s + K_d)/2;
    disp("Initial stationary eq.: r = " + r_eq + " ; K = " + K_eq)
else
    % ============== 2rd) Firms ========================================== %
    w = (1-alpha)*(alpha/(r_open_econ+dep))^(alpha/(1-alpha)); % equilibrium wages
    % ============== 3nd) Fiscal Policy ================================== %

    % 3.1 - Lump-sum transfers
    lump_sum = tau*w*Labor;   % Lump-Sum transfer to all households

    % ============== 4th) Households ===================================== %

    % Stationary eq.
    if specification.ind_fixed_effects == 0
        [g_final,con_final,~,~] = HJB(parameters,0,a_grid, ...
                                            borrowing_limit,r_open_econ,w,lump_sum,tau);
    else
        g_cell = cell(n_fe);
        con_cell = cell(n_fe);
        sav_cell = cell(n_fe);
        for i=1:n_fe
            [x,y,z] = HJB(parameters,omega_grid(i), ...
                                a_grid,borrowing_limit,r_open_econ,w,lump_sum,tau);
            g_cell{i} = x;
            con_cell{i} = y;
            sav_cell{i} = z;
        end
        g_final = 0;
        con_final = 0;
        sav_final = 0;
        for i = 1:n_fe
            g_final = g_final + omega_dist(i).*g_cell{i}.*domega_tilde(i);
            con_final = con_final + omega_dist(i).*con_cell{i}.*domega_tilde(i);
            sav_final = sav_final + omega_dist(i).*sav_cell{i}.*domega_tilde(i);
        end
    end
    r_eq = r_open_econ;
    a_marg_dist = g_final*dz_tilde; % marginal wealth dist
    K_eq = sum(a_grid.*a_marg_dist.*da_tilde);
end

%% ============== 5th) Equilibrium Objects ============================= %%

% 1) Aggregate Variables
Y = Labor*(alpha/(r_eq+dep))^(alpha/(1-alpha)); % output of the firm
C = (con_final.*g_final*dz_tilde)'*da_tilde; % aggregate consumption
S = Y-C; % aggregate savings
L = Labor; % Labor supply and demand
% Lump sum and capital have already been determined in the loop

% 2) Equilibrium prices 
w_eq = w;
% equilibrium/open econ interest rate determined in loop

% 3) Density 
% equilibrium density already determined in loop
a_marg_dist = g_final*dz_tilde; % marginal distribution of wealth

% 4) Store as a struct object
ss.output = Y;
ss.consumption = C;
ss.investment = S;
ss.labor = L;
ss.capital = K_eq;
ss.wage = w_eq;
ss.interest = r_eq;
ss.lump_sum = lump_sum;
ss.joint_dist = g_final;
ss.wealth_dist = a_marg_dist;

end




%=========================================================================
%                   FUNCTIONS
%=========================================================================


function [omega_dist,omega_grid,domega_tilde] = FE_grid(n_fe,sigma_omega)
    k = 0.7;

    % grid
    omega_grid = SymmetricPowerSpacedGrid(n_fe,k,2*sigma_omega,0);
    domega_tilde = grid_measure(omega_grid,n_fe);

    % Distribution
    pd = makedist("Normal","mu",0,"sigma",sigma_omega);
    % pd = makedist("GeneralizedPareto","k",200,"sigma",2,"theta",0);
    % pd = makedist("Gamma","a",1/2,"b",1);

    omega_dist = discrete_fe_dist(omega_grid,pd);

    omega_diag = spdiags(domega_tilde,0,n_fe,n_fe);     % I*n_z x I*n_z diagonal matrix
    omega_dist = omega_diag\omega_dist;
end


function dist = discrete_fe_dist(grid,pd)
    n = length(grid);
    dist = zeros(n,1);
    steps = diff(grid);
    for i=1:n
       if i == 1
            dist(1) = cdf(pd,grid(i) + steps(i)*1/2);
       elseif i == n
            dist(i) = 1 - cdf(pd,grid(i) - steps(i-1)*1/2);
       else
            dist(i) = cdf(pd,grid(i)+ steps(i)*1/2) - cdf(pd,grid(i)- steps(i-1)*1/2);
       end    
    end
end

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
end

