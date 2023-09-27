% ========================================================================
%                   Estimation
% ========================================================================
% Description: Employs a local or quasi-global approach to solve the SMM
%              minimization
%              
%
% Inputs:
%   • parameters:
%       - struct with all the parameters of the HA model
%   • settings:
%       - settings.algorithm (if 1, do global approach. If 2, local
%       search with Nelder-Mead Downhill simplex algorithm)
%       - settings.multishare (how many of the guesses to select
%       - settings.n_multi (number of guesses)
%       - settings.guess_min (lower bounds for guesses)
%       - settings.guess_max (upper bounds for guesses)
%   • Sample moments: data moments to compare
%
%
% Author: Rafael Lincoln/ Tomas Martinez (adapted from his Julia code)
% ========================================================================


function [par_sol, nfeval] = estimation(parameters,settings,specification,sample_moments,country,initial_guess)
    % This function performs parameter estimation using various optimization algorithms.
    % Check if the optional parameter is provided
    if nargin < 6
        initial_guess = 0;
    end
    
    % HOUSEKEEPING
    % Define the objective function using a lambda function handle
    f = @(x0) objective_fct(x0,parameters,specification,settings,sample_moments);

    % MULTI START WITH SIMPLEX
    if settings.algorithm == 1
        % Perform multistart optimization using the specified algorithm settings
        [multisol, allsol] = multiStart(f, settings);
        
        % Return the parameter solution and number of function evaluations
        par_sol = multisol.par_sol;
        nfeval = multisol.nfeval;

        % save solutions
        save("estimation_output\global_solution_estimation_" + country + ".mat","allsol","multisol")

    % NELDER MEAD SIMPLEX
    elseif settings.algorithm == 2
        % Perform optimization using Nelder-Mead simplex algorithm
        [opt_res,~,~,output] = fminsearch(f, initial_guess, optimset('PlotFcns','optimplotfval','TolX', 1e-3, 'TolFun', 1e-3, ...
            'MaxFunEvals', 5000, 'MaxIter', 2000, 'Display', 'off'));
        
        disp(output.message)

        % Return the parameter solution and number of function evaluations
        par_sol = opt_res;
        nfeval = output.iterations;

        % save locally
        save("estimation_output\local_solution_estimation_" + country + ".mat","opt_res","output")
    end
end
