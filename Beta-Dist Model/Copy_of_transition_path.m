%=========================================================================
%                   Transition path after a MIT Shock
%=========================================================================
% Description: Solves the transition path for prices and aggregate
% variables after a MIT-shock to the borrowing constraint
% 
%
%=========================================================================

function [Y_t,C_t,I_t,K_t,gg_t,a_marg_t,r_t,w_t,lump_sum_t] = Copy_of_transition_path(...
    a_grid,bl_path,permanent_income,transitory_income,closed_economy,ind_fixed_effects)

% global variables
global alpha n_p n_t dep tau n_fe

%=========================================================================
%    1)  Income Process and grids (time, asset, productivity, # of FE)
%=========================================================================

% # of grid points
n_z = n_p*n_t;

% Persistence of shocks
beta_p = permanent_income(1);
beta_t = transitory_income(1);

% Variance of shocks
sigma_p = permanent_income(2);
sigma_t = transitory_income(2);

% Arrival of shocks
lambda_p = permanent_income(3);
lambda_t = transitory_income(3);

% Income process CTMC matrix + grid
[lambda_z,z_grid] = inc_matrix(beta_p,beta_t,sigma_p,sigma_t,lambda_p,lambda_t);
dz_tilde = grid_measure(z_grid,n_z); % Define measure for the productivity grid for integration

% Bi-dimensional grids and measures
N_a = length(a_grid);         % length of asset grid
da_tilde = grid_measure(a_grid,N_a); % measure
% Measure for the AxZ grid
mu_grid = kron(dz_tilde,da_tilde);              % I*n_z x 1 (kronecker product)
grid_diag = spdiags(mu_grid,0,N_a*n_z,N_a*n_z); % I*n_z x I*n_z diagonal matrix


% sparse income matrix
Lambda_z = Build_Lambda_matrix(lambda_z,n_z,N_a);


% Time grid
N_t = 300; % Final period - by this time, the system should've converged to the new SS
T = 75;
dt = T/N_t;

% 1.3 - individual fixed effects
sigma_omega = 0.5;

% Obtain discrete distribution of FE, if applies
if ind_fixed_effects == 1
    % Discrete dist + grid + measure
    [omega_dist,omega_grid,domega_tilde] = FE_grid(n_fe,sigma_omega); 
end


% ===============================================================
% 1) Define boundary conditions: initial distribution and
%  terminal value function for the MFG system + intiial and final
%  values of capital for initial guess of the trajectory of
%  prices
% ===============================================================

% 1.1 - Initial Stationary Equilibrium => initial density
[~,C_ini,~,K_ini,L_ini,~,~,~,g_ini,a_marg_dist_ini] = stationary_equilibrium(permanent_income,transitory_income,a_grid, ...
    bl_path(1),closed_economy,ind_fixed_effects);

% 1.2 - Final stationary equilibrium
[~,~,~,K_fin,L_fin,lump_sum_fin,r_fin,w_fin,~,~] = stationary_equilibrium(permanent_income,transitory_income,a_grid, ...
    bl_path(N_t),closed_economy,ind_fixed_effects);

% Terminal condition is included inside the loop

% ===================================================================
% 2nd step - Guess path for capital in the transition path
% ===================================================================

L = (L_ini + L_fin)/2; % Labor = ergodic dist of z
K_t = linspace(K_ini,K_fin,N_t); % Intial guess for the path of capital
K_new=K_t;       %This is just preallocation. The values will not be used.

% Loop parameters
maxit = 1000;
convergence_criterion = 10^(-5);
relax= 0.5;

if ind_fixed_effects == 0
    % Terminal condition for the value function:
    [~,~,~,V_fin] = HJB(permanent_income,transitory_income,0,a_grid, ...
        bl_path(N_t),r_fin,w_fin,lump_sum_fin);

    % 3rd step - Solve time-dep. HJB given path of prices + terminal cond.

    % preallocation
    C_t = zeros(N_t,1); % aggregate consumption
    marginal_dists = cell(N_t,1);
    for it = 1:maxit
        fprintf('ITERATION = %d\n',it);
        % Guessing path for prices given path of capital
        r_t = alpha*(K_t.^(alpha-1)).*L.^(1-alpha) - dep;
        w_t = (1-alpha)*(K_t.^(alpha)).*L.^(-alpha);
        lump_sum_t = tau*w_t*L;
        % ============== Block 1: Solve backwards value functions =============
        [L_t,con_t] = vfi_time(r_t,w_t,lump_sum_t,0,V_fin,Lambda_z,a_grid,z_grid,N_t, ...
            bl_path,dt);

        % ================ Block 2: Solve distributions forward ===============
        gg = KF_time(g_ini,L_t,N_a,n_z,N_t,dt,grid_diag);

        % Take distributions and calculate objects of interest + update guess
        for n = 1:N_t
            if n == 1
                K_new(n) = K_ini; % Initial Capital
                C_t(1) = C_ini; % Initial Consumption
                marginal_dists{1} = a_marg_dist_ini; % Initial marginal dist
            else
                % reescale in order for integral to sum 1
                gg_reesc = grid_diag\gg(:,n);
                g = reshape(gg_reesc,N_a,n_z);
                % marginal wealth dist
                a_marg_dist = g*dz_tilde;
                % Aggregate objects of the economy
                marginal_dists{n} = a_marg_dist; % marginal dist
                K_new(n) = sum(a_grid.*a_marg_dist.*da_tilde); % aggregate capital
                C_t(n) = (con_t(:,:,n).*g*dz_tilde)'*da_tilde; % Aggregate Consumption
            end
        end

        % Check convergence and if not, update guess:
        fprintf('Maximum change in capital is %.8f\n',max(abs(K_t-K_new)));
        if max(abs(K_t-K_new))<convergence_criterion
            disp("Transition path converged!")
            break
        end

        % Update guess
        K_t=relax.*K_t+(1-relax).*K_new;
    end

else % CASE WITH INDIVIDUAL FIXED EFFECTS
    V_term = cell(n_fe);
    for i=1:n_fe % for each kind of individual calculate terminal condition for the value function
        [~,~,~,V_term{i}] = HJB(permanent_income,transitory_income,omega_grid(i),a_grid, ...
            bl_path(N_t),r_fin,w_fin,lump_sum_fin);
    end
    
    % 3rd step - Solve time-dep. HJB given path of prices + terminal cond. 
    
    % preallocation
    C_t = zeros(N_t,1); % aggregate consumption
    marginal_dists = cell(N_t,1);
    gg_t_ind = cell(n_fe,1);
    con_t_ind = cell(n_fe,1);
    for it = 1:maxit
        fprintf('ITERATION = %d\n',it);
        % Guessing path for prices given path of capital
        r_t = alpha*(K_t.^(alpha-1)).*L.^(1-alpha) - dep;
        w_t = (1-alpha)*(K_t.^(alpha)).*L.^(-alpha);
        lump_sum_t = tau*w_t*L;

        for i=1:n_fe
            % ============== Block 1: Solve backwards value functions =============
            [L_t,con_t_ind{i}] = vfi_time(r_t,w_t,lump_sum_t,omega_grid(i),V_term{i},Lambda_z,a_grid,z_grid,N_t, ...
            bl_path,dt);

            % ============== Block 2: Solve distributions forward ===============
            gg_t_ind{i} = KF_time(g_ini,L_t,N_a,n_z,N_t,dt,grid_diag);
        end

        g_final = 0;
        con_final = 0;
        % calculate aggregate consumption, capital and distribution
        for i = 1:n_fe
            g_final = g_final + omega_dist(i).*gg_t_ind{i}.*domega_tilde(i);
            con_final = con_final + omega_dist(i).*con_t_ind{i}.*domega_tilde(i);
        end


        % Take distributions and calculate objects of interest + update guess
        for n = 1:N_t
            if n == 1
                K_new(n) = K_ini; % Initial Capital
                C_t(1) = C_ini; % Initial Consumption
                marginal_dists{1} = a_marg_dist_ini; % Initial marginal dist
            else
                % reescale in order for integral to sum 1
                gg_reesc = grid_diag\g_final(:,n);
                g = reshape(gg_reesc,N_a,n_z);
                % marginal wealth dist
                a_marg_dist = g*dz_tilde;
                % Aggregate objects of the economy
                marginal_dists{n} = a_marg_dist; % marginal dist
                K_new(n) = sum(a_grid.*a_marg_dist.*da_tilde); % aggregate capital
                C_t(n) = (con_final(:,:,n).*g*dz_tilde)'*da_tilde; % Aggregate Consumption
            end
        end

        % Check convergence and if not, update guess:
        fprintf('Maximum change in capital is %.8f\n',max(abs(K_t-K_new)));
        if max(abs(K_t-K_new))<convergence_criterion
            break
        end

        % Update guess
        K_t=relax.*K_t+(1-relax).*K_new;
    end
    gg = g_final;
end

% ======= Other Aggregate Objects of the economy ======== %

% Output
Y_t = L.*(alpha./(r_t+dep)).^(alpha/(1-alpha));
% Investment/Savings
I_t = Y_t - C_t';
% distributions
gg_t= gg;
% marginal distributions
a_marg_t = marginal_dists;

% Net wealth to output Ratio
% wealth_output_t = K_t./Y_t;
% Debt-to-output

end




%=========================================================================
%                   FUNCTIONS
%=========================================================================


function [L_t,con_t] = vfi_time(r_t,w_t,lump_sum_t,fix_eff,V_T,Lambda_z,a_grid,z_grid,N_t, ...
                                bl_path,dt)
    global gamma

    N_a = length(a_grid);         % length of asset grid
    n_z = length(z_grid);          % length of productivity grid
    % path for value functions with terminal condition
    v = zeros(N_a,n_z,N_t); % value functions in time
    L_t = cell(N_t,1); % infinitesimal generators in time
    con_t = zeros(N_a,n_z,N_t); % consumption policy Functions

    % Preferences
    if gamma == 1
        u = @(x) log(x);
    else
        u = @(x) (x.^(1-gamma))./(1-gamma);
    end


    % v(:,:,N_t) = V_T; % set terminal condition
    V = V_T; % set terminal condition
    for n=N_t:-1:1
        % update value function
        v(:,:,n) = V;
        % index associated with the borrowing limit
        bl_n = bl_path(n);
        [~,idx_n] = min((a_grid-bl_n).^2);
        % build A Matrix
        [A,con] = Build_A_matrix(V,w_t(n),r_t(n),lump_sum_t(n),fix_eff,a_grid,idx_n,z_grid);
        u_n = u(con); % u at optimal c
        % Infinitesimal generator
        L_gen = A + Lambda_z;
        % Calculating value function next period
        V_next = update_V(V,u_n,L_gen,n_z,N_a,dt);
        % ========= Store objects of interest ========= %
        % infinitesimal generator
        L_t{n} = L_gen;
        % Consumption Policy function
        con_t(:,:,n) = con;
        % Savings policy function
        % sav = (1-tau)*w_t(n).*exp(zz + fix_eff) + r_t(n).*aa + lump_sum_t(n) - con;
        % Value function
        V = V_next;
    end
end

function gg = KF_time(g_ini,L_t,N_a,n_z,N_t,dt,grid_diag)
    % preallocation fo distributions in time
    gg = zeros(N_a*n_z,N_t); % distributions in time
    % Initial Condition
    gg0 = reshape(g_ini,N_a*n_z,1);
    % Reescale back to calculate correct the distribution
    gg_tilde_ini = grid_diag*gg0;
    gg(:,1) = gg_tilde_ini;
    for n = 1:N_t
        LT = L_t{n}'; % Infinitesimal Generator
        % Implicit method: g^{t+∆t} = (I - ∆t*L')^{-1}*(g^{t})
        gg(:,n+1) = (speye(N_a*n_z)-dt*LT)\gg(:,n);
    end    
end


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



% Enhancements to be done in the future

% Linear adjustment path
% bl_path = linspace(a_min_initial,a_min_final,N_t); 
% 
% % asset grid: merge normal asset grid + bl_path
% 
% a_grid = power_grid(borrowing_limit,amax,I);
% a_append = flip(bl_path)';
% a_grid = [a_append(1:N_t-1); a_grid];
% N_a = length(a_grid);



% Non-linear adjustment path
% nu = 0.2;
% bl_path = zeros(N,1);
% bl_path(1) = a_min_initial; bl_path(N) = a_min_final;
% for n=1:N-1
%     bl_path(n+1) = bl_path(n) + dt*nu*(a_min_final-bl_path(n));
% end


