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

settings.algorithm = 1;                                     % Algorithm chosen. if 1, global. If 2, local search
settings.multishare = 0.05;                                 % use x % of the guesses on 2nd stage
settings.n_multi = 1000;                                    % total number of guesses
% order: rho, delta_rho_1,delta_rho_2, tau, borrowing_limit
settings.guess_min = [0.01,1e-4,1e-4,0,-4.5];                    % lower boundaries for global search algorithm
settings.guess_max = [0.10,1e-2,1e-2,0.5,-1];                     % upper boundaries for global search algorithm
settings.weights = [1/4 1/4 1/4 1/10 0 1/20 1/20 1/20];  % Weights for the loss function (must sum up to 1)

% Notes on the chosen weights for the weighting matrix:
% - P20P40: weight is zero as this kind of model with unsecured debt makes
%   it impossible to match this quintile properly. (0)
% - P0P20: Assign higher weight as it this quintile concentrates the higher
%   share of wealthy HtM and inidividuals with higher MPCs, important for
%   the PE effect of the credit shock (double of the the other three
%   remaining) (2/5 x 1/4 = 1/10) 
% - P40P60, P60P80, P80P100: Weights are non-negative, as these moments are
%   important. However, less important than the first moment which has
%   double the weight (1/20)


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

% calibration = [0.0852,0.0003,0.1341,-3];
% borrowing_limit = calibration(4);

% ============================ Asset grid =============================== %

% Asset grid
I = 100;                            % # of points in the asset grid
a_grid = power_grid(borrowing_limit,parameters.amax,I);
n_add = 21;
bl_after = -5;
a_append = linspace(bl_after,borrowing_limit,n_add)';
a_grid = [a_append(1:n_add-1); a_grid];


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

% settings.algorithm = 2;
% [par_sol, nfeval] = estimation(parameters,settings,specification,sample_moments,country,calibration);

[par_sol, nfeval] = estimation(parameters,settings,specification,sample_moments,country);

% 
% calibration = [0.0546,0.0001,0.1677,-0.7463];
% 
% tic;
% loss = objective_fct(calibration,parameters,specification,settings,sample_moments);
% toc;
% 

calibration = [0.0810,0.005,0.005,0.2194,-1.6764];

calibration = allsol.par_sol(:,5);

ss = stationary_equilibrium(calibration,parameters,specification,a_grid);

sim_moments = Moments(ss,parameters,a_grid);

[loss,decomp] = loss_fct(sample_moments,sim_moments,settings.weights);






