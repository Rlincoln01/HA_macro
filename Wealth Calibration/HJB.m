%=========================================================================
%                   HJB 
%=========================================================================
%
% 
%=========================================================================

%=========================================================================
%                 1)  Income Process and asset grid
%=========================================================================

function [g,con,sav,V] = HJB(parameters,fixed_effect,a_grid, ...
                                            borrowing_limit,r,w,lump_sum,rho)

% check if rho is provided separately or not
if nargin < 8
    % if not provided, I use the value from the parameters
    rho = parameters.rho;
end

% options
display_iterations = 0;

% parameters and calibration
gamma = parameters.gamma;
tau = parameters.tau;

%=========================================================================
%                 1)  Income Process and asset grid
%=========================================================================

% # of grid points
n_p = parameters.n_p;
n_t = parameters.n_t;
% total # of points of the productivity grid
n_z = n_p*n_t;

% Income process CTMC matrix + grid
[lambda_z,z_grid] = inc_matrix(parameters);
% measure 
dz_tilde = grid_measure(z_grid,n_z);

% ergodic log-productivity distribution
% z_dist = stationary_dist(lambda_z,n_z,dz_tilde);

% fixed effect
omega = fixed_effect;

% ============== Asset grid ============= %

% index associated with the borrowing limit
[~,idx] = min((a_grid-borrowing_limit).^2);

% Number of points
n_a = length(a_grid);

% Build varying measure for the asset grid
da_tilde = grid_measure(a_grid,n_a);


% ============= Rest ==================== %

% 4.1) Determine joint grid of assets and income
zz = ones(n_a,1)*z_grid';                   % Matrix IxM of productivity
aa = a_grid*ones(1,n_z);                  % Matrix IxM of assets


% 1) Preferences
if gamma == 1
    u = @(x) log(x);
else
    u = @(x) (x.^(1-gamma))./(1-gamma);
end


% 4.2) Determine Sparse Lambda_z matrix
sparse_lambda_z = Build_Lambda_matrix(lambda_z,n_z,n_a);

% intial guess value function
v0 = u(max(exp(1)^(-10),lump_sum + r.*aa + (1-tau)*w.*exp(omega + zz)))./rho; 

% set v as initial guess
v = v0; % v^{0} that will be an input on the A matrix

% 4.4) Iteration Parameters
maxit = 100;        % Maximum number of iterations in the HJB loop
crit = 10^(-6);     % Criterion HJB loop
delta = 1000;       % Time Step 

for n=1:maxit
    % Update matrix v
    V = v; 

    % == 1st) Build the A^n matrix == % 
    [A,c] = Build_A_matrix(V,w,r,lump_sum,omega,a_grid,idx,z_grid,parameters);    
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

    V = update_V(V,u_n,L,n_z,n_a,delta,rho);    

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
% fix one value of the matrix, otherwise it has no inverse (its singular)
n = n_a*n_z;
% i_fix = (n+1)/2;
i_fix = 1;
row = [zeros(1,i_fix-1),1,zeros(1,n-i_fix)];
LT(i_fix,:) = row;

% Set the array of zeroes and also fix a value
b = zeros(n,1);
b(i_fix) = 1;

% stationary dist
g_tilde = LT\b;
g_sum = g_tilde'*ones(n,1);
g_tilde = g_tilde./g_sum;

gg = grid_diag\g_tilde;

g = reshape(gg,n_a,n_z);


end