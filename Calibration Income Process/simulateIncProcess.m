% ========================================================================
%                   Simulation
% ========================================================================
% Description: Simulates the quarterly income path for n_sim individuals or
% T years, and returns two matrices with yearly log income and yearly
% income in levels
%
%
% Author: Rafael Lincoln / Kaplan et al. (2018) (adapted from their
%         fortran code with changes)
% ========================================================================


function yannsim = simulateIncProcess(parameters,settings,draws)


% ==== Settings of the simulation ==== %

dt = 0.25;                      % time step in quarters
nsim = settings.n_sim;          % Number of income path simulations
n_years = settings.t_sim;       % Number of years simulated
Tsim = floor(n_years/dt) + 1;   % Number or quarters simulated
n_t = settings.n_t;
n_p = settings.n_p;



% ============1) Preparating the Markov Chains for simulation =========== %

% 1.1 Obtain the income matrix of the income processes and stationary dists
[~,~,TM_p,p_grid,TM_t,t_grid] = inc_matrix(parameters,settings);
p_dist = max(stationary_dist(TM_p,n_p),0);
t_dist = max(stationary_dist(TM_t,n_t),0);

% 1.2 - Discretize Transition Matrix
p_eye = eye(n_p); 
t_eye = eye(n_t);
t_trans = dt*TM_t + t_eye;
p_trans = dt*TM_p + p_eye;

% 1.3 - Define these objects as quarterly MC transition matrices
TM_p_qtr = p_trans;
TM_t_qtr = t_trans;

for it = 1:(floor(1/dt)-1)
    TM_p_qtr = TM_p_qtr*p_trans;
    TM_t_qtr = TM_t_qtr*t_trans; 
end

% ========== 2) Simulating income for n individuals in t qtrs =========== %


% 2.1 Define cumulative time vector (in quarters)
cumtvec = (1:Tsim)*dt;

% 2.2 Call random numbers for simulation
t_rand = draws.t_shocks;
p_rand = draws.p_shocks;

% 2.3 Preallocate income for individual and period
tsimI = zeros(nsim,Tsim); % transitory income (log,qtr)
psimI = zeros(nsim,Tsim); % persistent income (log,qtr)
ysim = zeros(nsim,Tsim); % sum of both (log,qtr)
ylevsim = zeros(nsim,Tsim); % sum of both (levels,qtr)
yannlevsim = zeros(nsim,Tsim); % sum of both (levels,yr)
yannsim = zeros(nsim,Tsim); % sum of both (log,yr)

for in=1:nsim
     % Simulate income path in dt increments
     % First shocks are drawn from the stationary distribution
     tsimI(in,1) = gendist(t_dist',1,1,t_rand(in,1));
     psimI(in,1) = gendist(p_dist',1,1,p_rand(in,1));
     ysim(in,1) = t_grid(tsimI(in,1)) + p_grid(psimI(in,1));

    for it = 2:Tsim
        % now, take shocks and simulate next period using transition matrix
        tsimI(in,it) = gendist(t_trans(tsimI(in,it-1),:),1,1,t_rand(in,it));
        psimI(in,it) = gendist(p_trans(psimI(in,it-1),:),1,1,p_rand(in,it));
        ysim(in,it) = t_grid(tsimI(in,it)) + p_grid(psimI(in,it));
    end

    %aggregate to annual income
    ylevsim(in,:) = exp(ysim(in,:));
    for it = 1:5
        yannlevsim(in,it) = sum(ylevsim(in,cumtvec>4.0*(it-1) & cumtvec<=4.0*it));
    end
    yannsim(in,:) = log(yannlevsim(in,:));
end

end
