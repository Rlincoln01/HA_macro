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
settings.n_multi = 10000;                      % total number of guesses
% order: borrowing limit, rho, delta_rho
settings.guess_min = [0,0.025,1e-4];           % lower boundaries for global search algorithm
settings.guess_max = [-10,0.09,1e-2];          % upper boundaries for global search algorithm


% ======================= Model Specification =========================== %

specification.ind_fixed_effects = 0;           % Compute individual fixed effects
specification.closed_economy = 1;              % Compute closed economy interest rate
specification.disc_rate_heterogeneity = 1;     % Equilibrium with discount rate heterogeneity


% ======================= Parametrization =============================== %

parameters.sigma_t = 1.74;                     % Variance of transitory shocks (from income process calibration)
parameters.sigma_p = 1.53;                     % Variance of persistent shocks (from income process calibration)
parameters.beta_t = 0.761;                     % Persistence of transitory shocks (from income process calibration)
parameters.beta_p = 0.009;                     % Persistence of permanent shocks (from income process calibration)
parameters.lambda_p = 0.007;                   % Arrival of transitory shocks (from income process calibration)
parameters.lambda_t = 0.08;                    % Arrival of permanent shocks (from income process calibration)
parameters.sigma_omega = 0.5;                  % Variance of the individual fixed effect
parameters.gamma = 2;                          % Relative risk aversion
parameters.alpha = 0.33;                       % Share of capital in production
parameters.zeta = 1/180;                       % arrival of death shock
parameters.dep = 0.05;                         % depreciation
parameters.tau = 0.15;                         % Tax rate
parameters.r_open_econ = 0.03;                 % Open-economy interest rate
parameters.amax = 600;                         % Maximum asset value
parameters.pgrid_width = 5.69800327154487/2;   % Persistent grid width
parameters.tgrid_width = 3.78503497352302/2;   % Transitory grid width
parameters.tgrid_par = 0.792311684654520;      % Level of non-linearity of the grid
parameters.pgrid_par = 0.641191531992662;      % Level of non-linearity of the grid
parameters.n_t = 3;                            % # of points in the transitory inc grid
parameters.n_p = 11;                           % # of points in the persistent inc grid
parameters.n_rho = 5;                          % # of points on the discount rate grid
parameters.delta_rho = 0.005;                  % distance between equally-space grid of discount rates
parameters.n_fe = 3;                           % # of points in the persistent inc grid


