%=========================================================================
%                   Transition path after a MIT Shock
%=========================================================================
% Description: Solves the transition path for prices and aggregate
% variables after a MIT-shock to the borrowing constraint
% 
%
%=========================================================================

global borrowing_limit amax I alpha gamma

%=========================================================================
%                 1)  Income Process and grids
%=========================================================================

% 1) Preferences
if gamma == 1
    u = @(x) log(x);
else
    u = @(x) (x.^(1-gamma))./(1-gamma);
end

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
nz = length(z_grid);

% Time grid
N_t = 100;
T = 75;
dt = T/N_t;
N1 = N; % Final period - by this time, the system should've converged to the new SS

% ====================== Borrowing Limit MIT shock ====================== %

a_min_initial = -borrowing_limit;
a_min_final = -2;
% Linear adjustment path
bl_path = linspace(a_min_initial,a_min_final,N_t); 


% Asset grids
a_grids = cell(N_t,1); % for preallocation
a_grids{1} = power_grid(a_min_initial,amax,I);
for i=1:N_t-1
    a_grids{i+1} = unique([bl_path(i+1);a_grids{i}]);
end

if ind_fixed_effects == 0
    % ===================================================================
    % 1st step - Compute stationary equilibria before and after the shock
    % ===================================================================
    
    % 1.1 - Initial Stationary Equilibrium => initial density
    [~,~,~,K_ini,L_ini,lump_sum_ini,r_ini,w_ini,g_ini,~] = stationary_equilibrium(permanent_income,transitory_income,a_grids{1}, ...
        closed_economy,ind_fixed_effects);
    % 1.2 - Final stationary equilibrium
    [~,~,~,K_fin,L_fin,lump_sum_fin,r_fin,w_fin,g_fin,~] = stationary_equilibrium(permanent_income,transitory_income,a_grids{N_t}, ...
                        closed_economy,ind_fixed_effects);
    % Terminal condition for the value function:
    [~,~,~,V_fin,~] = HJB_KF_implicit(permanent_income,transitory_income,0,a_grids{N_t}, ...
                                            r_end,w_end,lump_sum_fin,compute_MPC);

    % ===================================================================
    % 2nd step - Guess path for capital in the transition path
    % ===================================================================
    
    L = (L_ini + L_fin)/2; % Labor = ergodic dist of z
    K_t = linspace(K_ini,K_fin,N_t); % Intial guess for the path of capital
   
    % ====================================================================
    % 3rd step - Solve time-dep. HJB given path of prices + terminal cond.
    % ====================================================================
   
    % preallocation
    v = cell(N_t,1); % value functions in time
    gg = cell(N_t+1,1); % distributions in time
    L_t = cell(N_t,1); % infinitesimal generators in time

    maxit = 1000;
    convergence_criterion = 10^(-5);
    
    v{N_t} = V_fin;

    for it = 1:maxit
        fprintf('ITERATION = %d\n',it);
        % Guessing path for prices given path fo capital
        r_t = alpha*(K_t.^(alpha-1)).*L.^(alpha) - dep;
        w_t = (1-alpha)*(K_t.^(alpha)).*L.^(alpha-1);
        lump_sum_t = tau*w_t*L;
        for n=N_t:-1:1
            % update value function
            v{n} = V;
            % build A Matrix
            [A,con] = Build_A(V,w_t(n),r_t(n),lump_sum_t(n),0,a_grids{n},z_grid);
            u_n = u(con); % u at optimal c
            % sparse income matrix
            Lambda_z = Build_Lambda_matrix(lambda_z,n_z,length(a_grids{n}));
            % Infinitesimal generator
            L = A + Lambda_z;
            % Calculating value function next period
            V_next = update_V(V,u_n,L,n_z,n_a,dt);
            % ========= Store objects of intrerest ========= %
            % infinitesimal generator
            L_t{n} = L;
            % Savings policy function
            zz = ones(n_a,1)*z_grid';     % Matrix IxM of productivity
            aa = a_grids{n}*ones(1,n_z);      % Matrix IxM of assets
            sav = (1-tau)*w.*exp(zz + omega) + r.*aa + lump_sum - con;
        end
        % Now, Iterate distribution Forward using the KF equation
        
    end
else    
end





















% Enhancements to be done in the future

% Non-linear adjustment path
% nu = 0.2;
% bl_path = zeros(N,1);
% bl_path(1) = a_min_initial; bl_path(N) = a_min_final;
% for n=1:N-1
%     bl_path(n+1) = bl_path(n) + dt*nu*(a_min_final-bl_path(n));
% end
