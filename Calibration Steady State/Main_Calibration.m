% =========================================================================
% Calibration: File for setting the algorithm that calibrates the model
%              
%
%
% • Parametrization:
% • Model specification:
%
%
%
% Author: Rafael Lincoln
% =========================================================================

% =================== Settings of the calibration ======================= % 

settings.algorithm = 1;                        % Algorithm chosen. if 1, global. If 2, local search
settings.multishare = 0.05;                    % use x % of the guesses on 2nd stage
settings.n_multi = 1000;                      % total number of guesses
% order: rho, delta_rho, tau, borrowing_limit
settings.guess_min = [0.01,1e-4,0,-4.5];          % lower boundaries for global search algorithm
settings.guess_max = [0.10,1e-2,0.5,0];          % upper boundaries for global search algorithm

% ======================= Model Specification =========================== %

specification.ind_fixed_effects = 0;           % Compute individual fixed effects
specification.closed_economy = 1;              % Compute closed economy interest rate
specification.disc_rate_heterogeneity = 1;     % Equilibrium with discount rate heterogeneity


% ======================= Parametrization =============================== %

% other parameters - not used
parameters.sigma_omega = 0.5;                  % Variance of the individual fixed effect
parameters.r_open_econ = 0.03;                 % Open-economy interest rate
parameters.n_rho = 5;                          % # of points on the discount rate grid
parameters.n_fe = 3;                           % # of points in the FE grid
parameters.N_t = 300;                          % # of points in the time grid
parameters.T = 75;                             % termination period of the time grid    

% Income process calibration for Brazil
parameters.gamma = 1;                          % Relative risk aversion
parameters.alpha = 0.33;                       % Share of capital in production
parameters.dep = 0.025;                        % depreciation
parameters.r_open_econ = 0.03;                 % Open-economy interest rate
parameters.pgrid_width = 5.69800327154487/2;   % Persistent grid width
parameters.tgrid_width = 3.78503497352302/2;   % Transitory grid width
parameters.tgrid_par = 0.792311684654520;      % Level of non-linearity of the grid
parameters.pgrid_par = 0.641191531992662;      % Level of non-linearity of the grid
parameters.n_t = 3;                            % # of points in the transitory inc grid
parameters.n_p = 11;                           % # of points in the persistent inc grid
parameters.n_rho = 5;                          % # of points on the discount rate grid

% asset grid calibration
parameters.amin = -5;                          % minimum absolute value
parameters.amax = 600;                         % Maximum asset value
parameters.n_add = 21;                         % number of gridpoints to add after the borrowing limit

calibration = [0.0852,0.0003,0.1341,-3];
borrowing_limit = calibration(4);

% ============================ Asset grid =============================== %

% % Asset grid
% I = 100;                            % # of points in the asset grid
% a_grid = power_grid(borrowing_limit,parameters.amax,I);
% n_add = 21;
% bl_after = -5;
% a_append = linspace(bl_after,borrowing_limit,n_add)';
% a_grid = [a_append(1:n_add-1); a_grid];


% =========================== Calibration =============================== %

% select country
country = "BRA";

% sample moments of steady state
sample_moments_ss;
sample_moments = sample_moments.(country);

% Income process parameters
income_process_calibration;
income_parameters = inc_proc.(country);

moment_names = fieldnames(income_parameters);
n_f = length(moment_names);
F = zeros(n_f,1);

for i =1:n_f
    moment = moment_names{i};
    parameters.(moment) = income_parameters.(moment);
end

% Test with local search
% initial_guess = [0.05,0.005,0.15,-1.5];
% settings.algorithm = 2;
% [par_sol, nfeval] = estimation(parameters,settings,specification,sample_moments,country,initial_guess);
% 

[par_sol, nfeval] = estimation(parameters,settings,specification,sample_moments,country);


calibration = [0.0546,0.0001,0.1677,-0.7463];

tic;
loss = objective_fct(calibration,parameters,specification,settings,sample_moments);
toc;

ss = stationary_equilibrium(calibration,parameters,specification,a_grid);

sim_moments = Moments(ss,parameters,a_grid);
