% =========================================================================
% Main File:            
%
%
% • Parametrization:
% • Model specification:
%
%
%
% Author: Rafael Lincoln
% =========================================================================


% ======================= Model Specification =========================== %


ind_fixed_effects = 0;                   % Compute individual fixed effects
closed_economy = 1;                      % Compute closed economy interest rate
discount_rate_heterogeneity = 0;         % Equilibrium with discount rate heterogeneity
% compute_MPC % will remove and turn into a function


% ======================= Parametrization =============================== %

% Income Process
sigma_t = 1.74;                     % Variance of transitory shocks
sigma_p = 1.53;                     % Variance of persistent shocks
beta_t = 0.761;                     % Persistence of transitory shocks
beta_p = 0.009;                     % Persistence of permanent shocks
lambda_p = 0.007;                   % Arrival of transitory shocks
lambda_t = 0.08;                    % Arrival of permanent shocks
sigma_omega = 0.5;                  % Variance of the individual fixed effect

% Set vector of parameters for input of Stationary eq and transition path
permanent_income = [beta_p,sigma_p,lambda_p];
transitory_income = [beta_t,sigma_t,lambda_t];

% Preferences
global gamma rho
gamma = 2;                          % Relative risk aversion
rho = 0.05;                         % Discount rate

% Structural Parameters
global alpha zeta r_open_econ dep tau 
alpha = 0.33;                       % Share of capital in production
zeta = 1/180;                       % arrival of death shock
r_open_econ = 0.05227;              % open-economy interest rate
dep = 0.1;                          % depreciation
tau = 0.15;                         % Tax rate

% == Grid parametrization == %

global n_p n_t amax n_fe n_rho pgrid_width tgrid_width pgrid_par tgrid_par
% Permanent income grid
pgrid_width = 3.78503497352302;     % persistent grid width
pgrid_par = 0.641191531992662;      % Grid non-linearity parameter    
n_p = 11;                           % # of points in the persistent grid

% Transitory income grid
tgrid_width = 0.641191531992662;    % Transitory grid width
tgrid_par = 0.792311684654520;      % Grid non-linearity parameter    
n_t = 3;                            % # of points in the transitory grid

% Asset grid
borrowing_limit = 1;                % Borrowing limit for liquid wealth
amax = 600;                         % maximum asset holdings
I = 100;                            % # of points in the asset grid
a_grid = power_grid(-borrowing_limit,amax,I);


% Fixed effects grid
n_fe = 10;

% Discount rate grid
n_rho = 3;

% ======================= Call Steady State ============================= %

% Equilibrium distribution, output and prices
[a_marg_dist,Y,r_eq,w] = stationary_equilibrium(permanent_income,transitory_income, ...
                                               a_grid,closed_economy,ind_fixed_effects);

% Wealth statistics
[lorenz_df,wealth_shares,mean_liq_wealth] = wealth_stats(a_marg_dist,Y);

% ======================= Credit constraint MIT-shock =================== %





