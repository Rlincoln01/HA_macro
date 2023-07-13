%=========================================================================
%                   HJB Equation
%=========================================================================
% Description: Solves the system of IxJxK of the One-asset Aiyagari model 
% with a two-state income process bu using implicit method
% 
%
% Inputs: 
%
% - Income Calibration: 
%      inc_process = [permanent_params,trans_params,pgrid_width,
%                      tgrid_width]
% - Model Calibration
%      model_calibration = [asset_grid, model_params]
% - Fixed effect: \omega_{i} of the income process
%      
%
% where:
%   permanent_params = [n_p,beta_p,sigma_p,lambda_p]
%   trans_params = [n_t,beta_t,sigma_t,lambda_t]
%   asset_grid = [I, borrowing_limit]
%   model_params = [gamma, rho, alpha, dep, r, zeta, w, tau, lump_sum]
%
%
% Returns:
% - Stationary Distribution g
% - Policy functions c,s
% - Infinitesimal generator A, used to compute the stationary distribution
%   on the Kolmogorov Forward equation (?)
%=========================================================================

function [g,con,sav,V,MPCs] = HJB_KF_implicit(permanent_income,transitory_income,fixed_effect,a_grid, ...
                                            r,w,lump_sum,compute_MPC)

% Global variables defined in Main.m that are structural parameters of the model
global gamma rho zeta tau n_p n_t  

% options
display_iterations = 0;

if ~exist("compute_MPC","var")
    compute_MPC = 0;
end    

%=========================================================================
%                 1)  Income Process and asset grid
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
% measure 
dz_tilde = grid_measure(z_grid,n_z);

% ergodic log-productivity distribution
z_dist = stationary_dist(lambda_z,n_z,dz_tilde);

% Number of points
n_a = length(a_grid);

% Build varying measure for the asset grid
da_tilde = grid_measure(a_grid,n_a);

% fixed effect
omega = fixed_effect;


%=========================================================================
%                 2)  Model Calibration
%=========================================================================

% 1) Preferences
if gamma == 1
    u = @(x) log(x);
else
    u = @(x) (x.^(1-gamma))./(1-gamma);
end

% rest is determined by globals

%=========================================================================
%                 4)  HJB Iteration 
%=========================================================================

% 4.1) Determine joint grid of assets and income
zz = ones(n_a,1)*z_grid';                   % Matrix IxM of productivity
aa = a_grid*ones(1,n_z);                  % Matrix IxM of assets

% 4.2) Determine Sparse Lambda_z matrix

sparse_lambda_z = Build_Lambda_matrix(lambda_z,n_z,n_a);


% 4.3) Initial Guess <=> Terminal condition of the HJB
v0 = u(max(exp(1)^(-10),lump_sum + r.*aa + (1-tau)*w.*exp(omega + zz)))./(rho+zeta); 
% set v as initial guess
v = v0; % v^{0} that will be an input on the A matrix

% 4.4) Iteration Parameters
maxit = 100;        % Maximum number of iterations in the HJB loop
crit = 10^(-6);     % Criterion HJB loop
delta = 1000;       % Time Step 

% 4.5) Iteration Loop

for n=1:maxit
    % Update matrix v
    V = v; 

    % == 1st) Build the A^n matrix == % 
    [A,c] = Build_A(V,w,r,lump_sum,omega,a_grid,z_grid);
    % Max_{c} U(c)
    u_n = u(c);

    % == 2nd) Define the infinitesimal generator of the process
    L = A + sparse_lambda_z;

    % == 3rd) Check if it is a transition matrix (sum(row) = 1, for all rows)
    if max(abs(sum(L,2)))>10^(-9)
        disp('Improper Transition Matrix')
        break
    end    

    % == 4th) Update V == %

    V = update_V(V,u_n,L,n_z,n_a,delta);    

    % == 5th) Check for convergence == %

    Vchange = V - v;
    v = V;
   
    dist(n) = max(max(abs(Vchange)));
    if display_iterations == 1
        disp(['Value Function, Iteration ' int2str(n) ', max Vchange = ' num2str(dist(n))]);
    end
    if dist(n)<crit
        % disp('Value Function Converged, Iteration = ')
        % disp(n)
        break
    end    
end    
con = c;
sav = (1-tau)*w.*exp(zz + omega) + r.*aa + lump_sum - con;

%=========================================================================
%                       Kolmogorov-Forward Equation
%=========================================================================

% Obtain the auto-adjoint of the infinitesimal operator of V
LT = L';

% Measure for the AxZ grid
mu_grid = kron(dz_tilde,da_tilde);              % I*n_z x 1 (kronecker product)
grid_diag = spdiags(mu_grid,0,n_a*n_z,n_a*n_z); % I*n_z x I*n_z diagonal matrix

% ===================== For case where zeta = 0 ========================= %
% % fix one value of the matrix, otherwise it has no inverse (its singular)
% n = I*n_z;
% % i_fix = (n+1)/2;
% i_fix = 1;
% row = [zeros(1,i_fix-1),1,zeros(1,n-i_fix)];
% LT(i_fix,:) = row;
% 
% % Set the array of zeroes and also fix a value
% b = zeros(n,1);
% b(i_fix) = 1;
% 
% % stationary dist
% g_tilde = LT\b;
% g_sum = g_tilde'*ones(n,1);
% g_tilde = g_tilde./g_sum;
% 
% gg = grid_diag\g_tilde;
% 
% g = reshape(gg,I,n_z);

% ======================================================================= %

% Obtain the matrix of the death process:
B = -zeta.*speye(n_a*n_z);

% Discretize the distribution of newborns:
newborns = zeros(n_a*n_z,1);

% Newborns are born with approx zero assets and ergodic productivity distribution

[~,idx_a] = min(abs(a_grid)); % zero assets indexes

for m = 1:n_z   
    if m == 1
        index = idx_a;
        newborns(index) = z_dist(m);
    else    
        index = n_a*(m-1)+idx_a;
        newborns(index) = z_dist(m);
    end
    
end    

mass = da_tilde(idx_a);
newborns = (1/mass).*newborns; % dirac function * ergodic evalued at a = 0 on IxJxK vector

% ===========================================
% Solve the system 0 = L'*g + B*g + newborns
%               => g = -(L'+ B)^{-1}*newborns
%
% newborns = dirac*ergodic
% ===========================================
g_tilde = -(LT + B)\(zeta.*newborns);
g_sum = g_tilde'*ones(n_a*n_z,1); % reescale so it sums up to 1
g_tilde = g_tilde./g_sum;

gg = grid_diag\g_tilde;

g = reshape(gg,n_a,n_z);

%=========================================================================
%                       MPCs with Feynman-Kac formula
%=========================================================================

if compute_MPC == 1
    % Time grid
    qtr = 1; % # number of quarters of the MPC
    N = 100; % length of the time grid
    dt = qtr/N; % Step of the time grid
    
    % Terminal Condition
    Gamm_0 = zeros(n_z*n_a,1);
    
    % Re-estack consumption vector
    c_mpc = reshape(con,n_a*n_z,1);
    
    Gamm_before = Gamm_0;
    
    for i = 1:N
        % step iteration
        Gamm = ((1/dt)*speye(n_a*n_z) - L)\(c_mpc + (1/dt)*Gamm_before);
        % update for next step 
        Gamm_before = Gamm;
    end
    
    Gamm_sol = reshape(Gamm,n_a,n_z);
    
    % Compute MPCs by the derivatives
    
    MPCs = zeros(n_a-1,n_z);
    MPCs(1:n_a-1,:) = (Gamm_sol(2:n_a,:)-Gamm_sol(1:n_a-1,:))./daaf(1:n_a-1,:);
else
    MPCs = 0;
end

end

% a_adapted = a_grid(1:I-1);
% 
% 
% figure
% plot(a_adapted,MPCs(:,3)) % t = 0, p= -1,42
% hold on
% plot(a_adapted,MPCs(:,17))% t = 0, p= 0
% hold on
% plot(a_adapted,MPCs(:,31))% t = 0, p= 1,42
% legend(["t = 0, p= -1,42"," t = 0, p= 0","t = 0, p= 1,42"])
% xlim([-1.1,10]);
% 
% figure
% plot(a_adapted,MPCs(:,10)) % t = 0, p= -0,71
% hold on
% plot(a_adapted,MPCs(:,7))  % t = -0,94, p= 0
% legend(["t = 0, p= -0,71","t = -0,94, p= 0"])
% xlim([-1.1,10]);
% 
% zmin = z_grid(1);
% % MPCs 
% f = figure; 
% icut_a = 20;
% icut_z = 33;
% acut = a_grid(1:icut_a);
% zcut = z_grid(1:icut_z);
% concut = MPCs(1:icut_a,1:icut_z);
% set(gca,'FontSize',14)
% surf(acut,zcut,concut')
% view([45 25])
% xlabel('Wealth, $a$','FontSize',14,'interpreter','latex')
% ylabel('Productivity, $z$','FontSize',14,'interpreter','latex')
% zlabel('$MPC_{\omega,\tau}(a,z)$','FontSize',14,'interpreter','latex')
% xlim([amin max(acut)])
% ylim([zmin max(zcut)])
% 
