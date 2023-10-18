
%=========================================================================
%                   Transition path after a MIT Shock
%=========================================================================
% Description: Solves the transition path for prices and aggregate
% variables after a MIT-shock to the borrowing constraint
% 
%
%=========================================================================

function transition_allocation = transition_path_new(parameters,specification,a_grid,bl_path)

% Loop parameters
maxit = 1000;
convergence_criterion = 10^(-3);
relax= 0.5;

% Calibration
alpha = parameters.alpha;
dep = parameters.dep;
tau = parameters.tau;
rho = parameters.rho;

%=========================================================================
%       1)  Income Process and grids (time, asset and productivity)
%=========================================================================

% # of grid points
n_p = parameters.n_p;
n_t = parameters.n_t;
n_rho = parameters.n_rho;
delta_rho_1 = parameters.delta_rho_1;
delta_rho_2 = parameters.delta_rho_2;
n_z = n_p*n_t;

% Income process CTMC matrix + grid
[lambda_z,z_grid] = inc_matrix(parameters);
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
N_t = parameters.N_t; % Final period - by this time, the system should've converged to the new SS
T = parameters.T;
dt = T/N_t; % dt (time step) = quarter in the model


if specification.disc_rate_heterogeneity ==1
   % Calculate all discount rates
   rho_s = zeros(5,1);
   dist_max = delta_rho_1 + delta_rho_2;
   rho_s(1) = rho - dist_max;
   rho_s(5) = rho + dist_max;
   rho_s(2) = rho - delta_rho_1;
   rho_s(4) = rho + delta_rho_1;
   rho_s(3) = rho;

   % specification settings
   parameters_rho = parameters;
   specification_rho = specification;
   specification_rho.disc_rate_heterogeneity = 0;

   % Preallocation
   V_t_rhos = zeros(N_a,n_z,n_rho);
   K_new = zeros(N_t,n_rho);
   marginal_dists = zeros(N_a,N_t,n_rho);
   gg_t = zeros(N_a,n_z,N_t,n_rho);
   C_t = zeros(N_t,n_rho); % aggregate consumption
   Y_t = zeros(N_t,n_rho); % aggregate output
   I_t = zeros(N_t,n_rho); % aggregate investment
   K_t = zeros(N_t,n_rho); % aggregate capital
   r_t = zeros(N_t,n_rho);
   w_t = zeros(N_t,n_rho);
   lump_sum_t = zeros(N_t,n_rho);
   for j =1:n_rho
       fprintf("Calculating for individual %1.0f \n",j)
       rho_j = rho_s(j);
       parameters_rho.rho = rho_j;
       
       % Calculate initial and final steady state:
       
       % 1.1 - Initial Stationary Equilibrium => initial density
       iniSS = stationary_equilibrium(parameters_rho,specification_rho,a_grid, ...
           bl_path(1));

       % 1.2 - Final stationary equilibrium
       finSS = stationary_equilibrium(parameters_rho,specification_rho,a_grid, ...
           bl_path(N_t));

       % 1.3 - Terminal Condition
       [~,~,~,V_t_rhos(:,:,j)] = HJB(parameters_rho,0,a_grid, ...
           bl_path(N_t),finSS.interest,finSS.wage,finSS.lump_sum,rho_j);

       % Guess path for prices
       L = (iniSS.labor + finSS.labor)/2; % Labor = ergodic dist of z
       K_t(:,j) = linspace(iniSS.capital,finSS.capital,N_t); % Intial guess for the path of capital
       K_new(:,j)=K_t(:,j);       %This is just preallocation. The values will not be used.

       for it = 1:maxit
           fprintf('ITERATION = %d\n',it);
           % Guessing path for prices given path of capital
           r_t(:,j) = alpha*(K_t(:,j).^(alpha-1)).*L.^(1-alpha) - dep;
           w_t(:,j) = (1-alpha)*(K_t(:,j).^(alpha)).*L.^(-alpha);
           lump_sum_t(:,j) = tau*w_t(:,j)*L;
           % Calculate transition path for individual j:
           [L_t,con_t] = vfi_time(r_t(:,j),w_t(:,j),lump_sum_t(:,j), ...
               0,V_t_rhos(:,:,j),Lambda_z,a_grid,z_grid, ...
               N_t, ...
               bl_path,dt,parameters,rho_j);
           gg = KF_time(iniSS.joint_dist,L_t,N_a,n_z,N_t,dt,grid_diag);
           for n = 1:N_t
               if n == 1
                   K_new(n,j) = iniSS.capital; % Initial Capital
                   C_t(1,j) = iniSS.consumption; % Initial Consumption
                   marginal_dists(:,1,j) = iniSS.wealth_dist; % Initial marginal dist
                   gg_t(:,:,n,j) = iniSS.joint_dist;
               else
                   % reescale in order for integral to sum 1
                   gg_reesc = grid_diag\gg(:,:,n);
                   g = reshape(gg_reesc,N_a,n_z);
                   gg_t(:,:,n,j) = g;
                   % marginal wealth dist
                   a_marg_dist = g*dz_tilde;
                   % Aggregate objects of the economy
                   marginal_dists(:,n,j) = a_marg_dist; % marginal dist
                   K_new(n,j) = sum(a_grid.*a_marg_dist.*da_tilde); % Aggregate capital
                   C_t(n,j) = (con_t(:,:,n).*g*dz_tilde)'*da_tilde; % Aggregate Consumption
               end
           end


           % Check convergence and if not, update guess:
           fprintf('Maximum change in capital is %.8f\n',max(abs(K_t(:,j)-K_new(:,j))));
           if max(abs(K_t(:,j)-K_new(:,j)))<convergence_criterion
               fprintf("Transition Path for individual %1.0f converged!",j)
               break
           end

           % time_grid = (1:N_t)*dt;
           % plot(time_grid,r_t)
           % hold on
           % drawnow;
           % Update guess
           K_t(:,j)=relax.*K_t(:,j)+(1-relax).*K_new(:,j);
       end
       Y_t(:,j) = L.*(alpha./(r_t(:,j)+dep)).^(alpha/(1-alpha)); % Output
       I_t(:,j) = Y_t(:,j) - C_t(:,j); % Investment
   end
   K_t = mean(K_t,2);
   C_t = mean(C_t,2);
   r_t = mean(r_t,2);
   w_t = mean(w_t,2);
   I_t = mean(I_t,2);
   Y_t = mean(Y_t,2);
   lump_sum_t = mean(lump_sum_t,2);
   marginal_dists = mean(marginal_dists,3);
   a_marg_t = marginal_dists;
   gg_t = mean(gg_t,4);
else
    % 1.1 - Initial Stationary Equilibrium => initial density
    iniSS = stationary_equilibrium(parameters,specification,a_grid, ...
        bl_path(1));

    % 1.2 - Final stationary equilibrium
    finSS = stationary_equilibrium(parameters,specification,a_grid, ...
        bl_path(N_t));

    % Value function
    [~,~,~,V_T] = HJB(parameters,0,a_grid, ...
        bl_path(N_t),finSS.interest,finSS.wage,finSS.lump_sum);

    L = (iniSS.labor + finSS.labor)/2; % Labor = ergodic dist of z
    K_t = linspace(iniSS.capital,finSS.capital,N_t); % Intial guess for the path of capital
    K_new=K_t;       %This is just preallocation. The values will not be used.

    gg_t = zeros(N_a,n_z,N_t);
    marginal_dists = zeros(N_a,N_t);
    C_t = zeros(N_t,1);
    for it = 1:maxit
        fprintf('ITERATION = %d\n',it);
        % Guessing path for prices given path of capital
        r_t = alpha*(K_t.^(alpha-1)).*L.^(1-alpha) - dep;
        w_t = (1-alpha)*(K_t.^(alpha)).*L.^(-alpha);
        lump_sum_t = tau*w_t*L;
        % ============== Block 1: Solve backwards value functions =============
        [L_t,con_t] = vfi_time(r_t,w_t,lump_sum_t,0,V_T,Lambda_z,a_grid,z_grid,N_t, ...
            bl_path,dt,parameters);

        % ================ Block 2: Solve distributions forward ===============

        gg = KF_time(iniSS.joint_dist,L_t,N_a,n_z,N_t,dt,grid_diag);

        % Take distributions and calculate objects of interest + update guess
        for n = 1:N_t
            if n == 1
                K_new(n) = iniSS.capital; % Initial Capital
                C_t(1) = iniSS.consumption; % Initial Consumption
                marginal_dists(:,n) = iniSS.wealth_dist; % Initial marginal dist
                gg_t(:,:,n) = iniSS.joint_dist; % initial joint distribution
            else
                % reescale in order for integral to sum 1
                gg_reesc = grid_diag\gg(:,:,n);
                g = reshape(gg_reesc,N_a,n_z);
                gg_t(:,:,n) = g;
                % marginal wealth dist
                a_marg_dist = g*dz_tilde;
                % Aggregate objects of the economy
                marginal_dists(:,n) = a_marg_dist; % marginal dist
                K_new(n) = sum(a_grid.*a_marg_dist.*da_tilde); % aggregate capital
                C_t(n) = (con_t(:,:,n).*g*dz_tilde)'*da_tilde; % Aggregate Consumption
            end
        end

        % Check convergence and if not, update guess:
        fprintf('Maximum change in capital is %.8f\n',max(abs(K_t-K_new)));
        if max(abs(K_t-K_new))<convergence_criterion
            disp("Transition Path converged!")
            break
        end

        % time_grid = (1:N_t)*dt;
        % plot(time_grid,r_t)
        % hold on
        % drawnow;
        % Update guess
        K_t=relax.*K_t+(1-relax).*K_new;
    end
    % Output
    Y_t = L.*(alpha./(r_t+dep)).^(alpha/(1-alpha));
    % Investment/Savings
    I_t = Y_t - C_t';
    % marginal distributions
    a_marg_t = marginal_dists;
end

% ======= Other Aggregate Objects of the economy ======== %


% store transition paths as a struct
transition_allocation.output = Y_t;
transition_allocation.consumption = C_t';
transition_allocation.investment = I_t;
transition_allocation.capital = K_t;
transition_allocation.joint_dist = gg_t;
transition_allocation.wealth_dist = a_marg_t;
transition_allocation.interest = r_t;
transition_allocation.wage = w_t;
transition_allocation.lump_sum = lump_sum_t;

end

%=========================================================================
%                   FUNCTIONS
%=========================================================================


function [L_t,con_t] = vfi_time(r_t,w_t,lump_sum_t,fix_eff,V_T,Lambda_z,a_grid,z_grid,N_t, ...
                                bl_path,dt,parameters,rho)
    
    % check if rho is provided separately or not
    if nargin < 13
        % if not provided, I use the value from the parameters
        rho = parameters.rho;
    end

    gamma = parameters.gamma;       % relative risk aversion
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
        [A,con] = Build_A_matrix(V,w_t(n),r_t(n),lump_sum_t(n),fix_eff,a_grid,idx_n,z_grid,parameters);
        u_n = u(con); % u at optimal c
        % Infinitesimal generator
        L_gen = A + Lambda_z;
        % Calculating value function next period
        V_next = update_V(V,u_n,L_gen,n_z,N_a,dt,rho);
        % ========= Store objects of intrerest ========= %
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
% 

function gg = KF_time(g_ini,L_t,N_a,n_z,N_t,dt,grid_diag)
    % preallocation fo distributions in time
    gg = zeros(N_a*n_z,1,N_t+1); % distributions in time
    % Initial Condition
    gg0 = reshape(g_ini,N_a*n_z,1);
    % Reescale back to calculate correct the distribution
    gg_tilde_ini = grid_diag*gg0;
    gg(:,:,1) = gg_tilde_ini;
    for n = 1:N_t
        LT = L_t{n}'; % Infinitesimal Generator
        % Implicit method: g^{t+∆t} = (I - ∆t*L')^{-1}*(g^{t})
        gg(:,:,n+1) = (speye(N_a*n_z)-dt*LT)\gg(:,:,n);
    end    
end
