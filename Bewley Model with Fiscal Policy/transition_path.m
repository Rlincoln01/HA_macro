%=========================================================================
%                   Solve for Transition Path
%=========================================================================

% Loop parameters
maxit = 1000;
convergence_criterion = 10^(-5);
relax= 0.5;

% path for taxes
tau_0 = parameters.tau;
tau_fin = 0.1;
tau_path = linspace(tau_0,tau_fin,n_t);

% === Block 1: Equilibrium conditions + initial guess for prices === %

% Initial stationary equilibrium
[B_ini,r_ini] = stationary_equilibrium(parameters,ydist,tau_path(1));

% Final stationary equilibrium
[B_fin,r_fin] = stationary_equilibrium(parameters,ydist,tau_path(n_t));

% Terminal condition for the value function
V_fin = HJB(parameters,r,tau_path(n_t));

% initial guess for prices
r_t = linspace(r_ini,r_fin,n_t);


% preallocation
v = zeros(n_a,n_y,n_t);     % Prices
B_t = zeros(n_t,1);         % Bonds

V(:,:,n_t) = V_fin;

for it= 1:maxit
    fprintf('ITERATION = %d\n',it);
    % new guess based on changes
    r_t = tau_path.*(y*ydist)./B_t; % government budget constraint
    for n=n_t:-1:1
        v(:,:,n) = V;

end