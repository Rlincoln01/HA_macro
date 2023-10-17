%=========================================================================
%                   IRF decomposition
%=========================================================================
% Description: Calculates the IRF decomposition of an MIT shock into
% general and partial equilibrium effects
% 
%
% Inputs: 
%
% Outputs: 
%
%=========================================================================

function C_t_PE = PE_decomposition(parameters,specification,gamma_t,a_grid,iniSS)


%=========================================================================
%       0) Presetting grids and calibration
%=========================================================================


% Calibration
rho = parameters.rho;
delta_rho_1 = parameters.delta_rho_1;
delta_rho_2 = parameters.delta_rho_2;

% Income and asset grid
n_p = parameters.n_p;
n_t = parameters.n_t;
n_rho = parameters.n_rho;
n_z = n_p*n_t;
N_a = length(a_grid);         % length of asset grid
da_tilde = grid_measure(a_grid,N_a); % measure

% Income process CTMC matrix + grid
[lambda_z,z_grid] = inc_matrix(parameters);
dz_tilde = grid_measure(z_grid,n_z); % Define measure for the productivity grid for integration
mu_grid = kron(dz_tilde,da_tilde);              % I*n_z x 1 (kronecker product)
grid_diag = spdiags(mu_grid,0,N_a*n_z,N_a*n_z); % I*n_z x I*n_z diagonal matrix

% sparse income matrix
Lambda_z = Build_Lambda_matrix(lambda_z,n_z,N_a);

% Time grid
N_t = parameters.N_t; % Final period - by this time, the system should've converged to the new SS
T = parameters.T;
dt = T/N_t; % dt (time step) = quarter in the model

% If there is discount rate heterogeneity, then calculate all discount rates
if specification.disc_rate_heterogeneity ==1
    rho_s = zeros(5,1);
    dist_max = delta_rho_1 + delta_rho_2;
    rho_s(1) = rho - dist_max;
    rho_s(5) = rho + dist_max;
    rho_s(2) = rho - delta_rho_1;
    rho_s(4) = rho + delta_rho_1;
    rho_s(3) = rho;
end


%=========================================================================
%       Compute Partial Equilibrium effects
%=========================================================================

% 1.1 - Take as inputs the path Γ_t = {θ^i_t,θ^{-i}_t}
bl_path = gamma_t.bl_t;
w_t = gamma_t.wages_t;
r_t = gamma_t.interest_t;
lump_sum_t = gamma_t.lump_sum_t;

% 1.2 - Compute terminal condition(s) under partial equilibrium
if specification.disc_rate_heterogeneity == 1
    V_t_rhos = zeros(N_a,n_z,n_rho);
    for j =1:n_rho
        rho_j = rho_s(j);
        [~,~,~,V_t_rhos(:,:,j)] = HJB(parameters,0,a_grid, ...
        bl_path(N_t),r_t(N_t),w_t(N_t),lump_sum_t(N_t),rho_j);
    end
else
    [~,~,~,V_T] = HJB(parameters,0,a_grid, ...
        bl_path(N_t),r_t(N_t),w_t(N_t),lump_sum_t(N_t));
end


% 1.3 Solve time-dependent HJB and KF equations for PE

disp("Calculating Partial Equilibrium effects for the path Γ_t...")
if specification.disc_rate_heterogeneity == 1
gg_t = zeros(N_a,n_z,N_t,n_rho); % pre-allocation
C_t_PE = zeros(N_t,n_rho); % aggregate consumption
    for j=1:n_rho
        fprintf("Calculating for individual %1.0f \n",j)
        rho_j = rho_s(j);
        [L_t,con_t] = vfi_time(r_t,w_t,lump_sum_t,0,V_t_rhos(:,:,j),Lambda_z,a_grid,z_grid,N_t, ...
            bl_path,dt,parameters,rho_j);
        gg = KF_time(iniSS.joint_dist,L_t,N_a,n_z,N_t,dt,grid_diag);
        for n = 1:N_t
            if n == 1
                C_t_PE(1,j) = iniSS.consumption; % Initial Consumption
                gg_t(:,:,n,j) = iniSS.joint_dist;
            else
                % reescale in order for integral to sum 1
                gg_reesc = grid_diag\gg(:,:,n);
                g = reshape(gg_reesc,N_a,n_z);
                gg_t(:,:,n,j) = g;
                % Aggregate objects of the economy
                C_t_PE(n,j) = (con_t(:,:,n).*g*dz_tilde)'*da_tilde; % Aggregate Consumption
            end
        end
    end
    C_t_PE = mean(C_t_PE,2); % aggregate consumption
else
    gg_t = zeros(N_a,n_z,N_t);
    C_t_PE = zeros(N_t,1);

    % backward solve value functions
    [L_t,con_t] = vfi_time(r_t,w_t,lump_sum_t,0,V_T,Lambda_z,a_grid,z_grid,N_t, ...
        bl_path,dt,parameters);

    % solve forward distributions
    gg = KF_time(iniSS.joint_dist,L_t,N_a,n_z,N_t,dt,grid_diag);

    for n = 1:N_t
        if n == 1
            C_t_PE(n) = iniSS.consumption; % Initial Consumption
            gg_t(:,:,n) = iniSS.joint_dist;
        else
            % reescale in order for integral to sum 1
            gg_reesc = grid_diag\gg(:,:,n);
            g = reshape(gg_reesc,N_a,n_z);
            gg_t(:,:,n) = g;
            % Aggregate objects of the economy
            C_t_PE(n) = (con_t(:,:,n).*g*dz_tilde)'*da_tilde; % Aggregate Consumption
        end
    end
end    
disp("Done!")

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
