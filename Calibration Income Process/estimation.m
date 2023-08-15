% ========================================================================
%                   Estimation
% ========================================================================
% Description: Employs a local or quasi-global approach to solve the SMM
%              minimization
%              
%
% Inputs:
%   • settings:
%       - settings.algorithm (if 1, do global approach. If 2, local
%       search with Nelder-Mead Downhill simplex algorithm)
%       - settings.multishare (how many of the guesses to select
%       - settings.n_multi (number of guesses)
%       - settings.guess_min (lower bounds for guesses)
%       - settings.guess_max (upper bounds for guesses)
%   • draws:
%       - draws for transitory and income processes (use seed for
%       comparability)
%   • Sample moments: data moments to compare
%
%
% Author: Rafael Lincoln/ Tomas Martinez (adapted from his Julia code)
% ========================================================================


function [par_sol, nfeval] = estimation(settings, draws, sample_moments,initial_guess)
    % This function performs parameter estimation using various optimization algorithms.
    % Check if the optional parameter is provided
    if nargin < 4
        initial_guess = 0;
    end
    
    % HOUSEKEEPING
    % Define the objective function using a lambda function handle
    f = @(x0) objective_fct(x0,settings,draws,sample_moments);

    % MULTI START WITH SIMPLEX
    if settings.algorithm == 1
        % Perform multistart optimization using the specified algorithm settings
        [multisol, allsol] = multiStart(f, settings);
        
        % Return the parameter solution and number of function evaluations
        par_sol = multisol.par_sol;
        nfeval = multisol.nfeval;

    % NELDER MEAD SIMPLEX
    elseif settings.algorithm == 2
        % Perform optimization using Nelder-Mead simplex algorithm
        [opt_res,~,~,output] = fminsearch(f, initial_guess, optimset('PlotFcns','optimplotfval','TolX', 1e-6, 'TolFun', 1e-3, ...
            'MaxFunEvals', 5000, 'MaxIter', 2000, 'Display', 'off'));
        
        disp(output.message)

        % Return the parameter solution and number of function evaluations
        par_sol = opt_res;
        nfeval = output.iterations;
    end
end
