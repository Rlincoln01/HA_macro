% ========================================================================
%                   Multistart estimation
% ========================================================================
% Description: This function performs multi-start optimization with 
%              two stages.
%
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
%
%
% Author: Rafael Lincoln/ Tomas Martinez (adapted from his Julia code)
% ========================================================================


function [bestsol, allsol] = multiStart(f, algSettings)

    guess_min = algSettings.guess_min;
    guess_max = algSettings.guess_max;

    select = ceil(algSettings.multishare * algSettings.n_multi); % select % of all multistart
   
    % Generate random guesses using Sobol sequences
    np = length(guess_min);
    s = sobolset(np, 'Skip', 1e3, 'Leap', 1e2); % Initialize Sobol sequence with skip and leap values
    s = scramble(s, 'MatousekAffineOwen'); % Scramble the sequence
    % scale the points within the parameter bounds
    guesses = net(s,algSettings.n_multi).*(guess_max-guess_min) + guess_min; 

    % First stage of the multistart algorithm
    fprintf('Starting first stage of Multistart\n');
    start_time = tic;
    objec_val = zeros(1, algSettings.n_multi);
    WaitMessage = parfor_wait(algSettings.n_multi, 'Waitbar', true);
    parfor j = 1:algSettings.n_multi
        WaitMessage.Send;
        objec_val(j) = f(guesses(j,:)); 
    end
    %Destroy the object. 
    WaitMessage.Destroy
    timeM = toc(start_time);

    [~, index] = sort(objec_val); % Select index of best guesses based on objective
    opt_guess = guesses(index(1:select),:); % Select optimal guesses
    fprintf('First stage done in %1.2f minutes.\n', timeM / 60);

    % Save output of the first stage
    save("estimation_output\first_stage_global_est.mat","guesses","objec_val","opt_guesses")

    fprintf('\nStarting second stage of Multistart\n');
    start_time = tic;

    fvals = zeros(1, select);
    par_sol = zeros(np, select);
    nfeval = zeros(1, select);

    parfor i = 1:select
        fprintf('Multistart iteration: %d\n', i);
        [opt_res,fval,~,output] = fminsearch(f, opt_guess(i,:), optimset('TolX', 1e-6, 'TolFun', 1e-3, ...
            'MaxFunEvals', 500, 'MaxIter', 200, 'Display', 'off'));
        nfeval(i) = output.iterations;
        par_sol(:, i) = opt_res;
        fvals(i) = fval;
    end
    timeM = toc(start_time);
    fprintf('Second stage done in %1.2f minutes.\n', timeM / 60);

    allsol = struct('fvals', fvals, 'par_sol', par_sol, 'nfeval', nfeval);
    
    bestsol = struct('fvals', fvals(1), 'par_sol', par_sol(:, 1), ...
        'nfeval', sum(nfeval)+ algSettings.n_multi);
end
