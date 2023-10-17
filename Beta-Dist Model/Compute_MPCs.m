%=========================================================================
%                       MPCs with Feynman-Kac formula
%=========================================================================
%
%
% 
%=========================================================================

function [MPCs_windfall] = Compute_MPCs(qtr,transfer,parameters,specification,fixed_effect,a_grid,...
                        borrowing_limit,r,w,lump_sum)


% calibration/parameters
rho = parameters.rho;
tau = parameters.tau;
n_p = parameters.n_p;
n_t = parameters.n_t;
n_rho = parameters.n_rho;
delta_rho_1 = parameters.delta_rho_1;
delta_rho_2 = parameters.delta_rho_2;
gamma = parameters.gamma;

% options
display_iterations = 0;

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
%                 1)  Income Process and asset grid
%=========================================================================

% # of grid points
n_z = n_p*n_t;

% Income process CTMC matrix + grid
[lambda_z,z_grid] = inc_matrix(parameters);

% fixed effect
omega = fixed_effect;

% ============== Asset grid ============= %

% index associated with the borrowing limit
[~,idx] = min((a_grid-borrowing_limit).^2);

% Number of points
n_a = length(a_grid);


% ============= Rest ==================== %

% 4.1) Determine joint grid of assets and income
zz = ones(n_a,1)*z_grid';                   % Matrix IxM of productivity
aa = a_grid*ones(1,n_z);                  % Matrix IxM of assets

% Build varying measure for the asset grid
[~, daf, ~] = grid_measure(a_grid,n_a);

daaf = daf*ones(1,n_z);           % Matrix IXJ of F-measures

% 1) Preferences
if gamma == 1
    u = @(x) log(x);
else
    u = @(x) (x.^(1-gamma))./(1-gamma);
end


% 4.2) Determine Sparse Lambda_z matrix
sparse_lambda_z = Build_Lambda_matrix(lambda_z,n_z,n_a);

% 4.4) Iteration Parameters
maxit = 100;        % Maximum number of iterations in the HJB loop
crit = 10^(-6);     % Criterion HJB loop
delta = 1000;       % Time Step 


if specification.disc_rate_heterogeneity == 1
MPCs_windfall = zeros(n_a,n_z,n_rho);        
    for j=1:n_rho
        disp("Solving Feynman-Kac equation for individual: " + j)
        % preallocation
        rho_j = rho_s(j);
        % intial guess value function
        v0 = u(max(exp(1)^(-10),lump_sum + r.*aa + (1-tau)*w.*exp(omega + zz)))./rho_j;
        % set v as initial guess
        v = v0; % v^{0} that will be an input on the A matrix
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
            V = update_V(V,u_n,L,n_z,n_a,delta,rho_j);
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
        % ======== Instantaneous MPC over a period tau ======== %
        % Time grid
        % qtr = 1; % # number of quarters of the MPC
        N = 300; % length of the time grid
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
        % ======== Compute the MPC of a windfall gain ======== %
        MPCs_windfall(:,:,j) = (interp2(zz,aa,Gamm_sol,zz,aa+transfer)-Gamm_sol)./transfer;
        MPCs_windfall(n_a,:,j) = MPCs_windfall(n_a-1,:,j);
    end
    MPCs_windfall = mean(MPCs_windfall,3);
else
    % intial guess value function
    v0 = u(max(exp(1)^(-10),lump_sum + r.*aa + (1-tau)*w.*exp(omega + zz)))./rho;
    % set v as initial guess
    v = v0; % v^{0} that will be an input on the A matrix
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
    % ======== Instantaneous MPC over a period tau ======== %
    % Time grid
    % qtr = 1; % # number of quarters of the MPC
    N = 300; % length of the time grid
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
    % ======== Compute the MPC of a windfall gain ======== %
    MPCs_windfall = (interp2(zz,aa,Gamm_sol,zz,aa+transfer)-Gamm_sol)./transfer;
    MPCs_windfall(n_a,:) = MPCs_windfall(n_a-1,:);
end    

end


    % Compute MPCs by the derivatives

    % MPCs = zeros(n_a,n_z);
    % MPCs(1:n_a-1,:) = (Gamm_sol(2:n_a,:)-Gamm_sol(1:n_a-1,:))./daaf(1:n_a-1,:);
    %
    % % approximate the MPCs for the last value
    % MPCs(n_a,:) = MPCs(n_a-1,:);
