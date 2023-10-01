% ========================================================================
%                   Objective function
% ========================================================================
% Description: Solves for the Steady State of the model, calculates the 
%              relevant moments and computes the loss function
%
% Inputs:
%   
%
%
% Author: Rafael Lincoln
% ========================================================================


function loss = objective_fct(calibration,parameters,specification,settings,sample_moments)

% Penalty function => fminsearch has no bound restrictions
penalty = penalty_fct(calibration,settings); %function checks whether param bounds have been violated

% solve for the SS and compute moments
if penalty == 0
    % Call asset grid based on borrowing limit
    borrowing_limit = calibration(4);
    a_grid = asset_grid(borrowing_limit,parameters);
    % Compute Steady State Equilibrium
    ss = stationary_equilibrium(calibration,parameters,specification,a_grid);

    % retrieve equilibrium objects for sanity checks
    K = ss.capital;
    Y = ss.output;
    C = ss.consumption;
    r_eq = ss.interest;

    sim_moments = Moments(ss,parameters,a_grid);

    if (r_eq <= 0) || (C>Y) || (K <= 0) 
        loss = 10e6;
    else 
        loss = loss_fct(sample_moments,sim_moments);
    end
else
    loss = 10e6;
end

end


% Penalty function

function penalty = penalty_fct(calibration,settings)
    % Preallocation
    penalty= 0;
    guess_min = settings.guess_min;
    guess_max = settings.guess_max;
    % If a parameter is below the minimum values
    condition_min = calibration >= guess_min;
    penalty = penalty + 10000*(min(min(condition_min),0)).^2;
    % If a parameter is above the maximum values
    condition_max = calibration <= guess_max;
    penalty = penalty + 10000*(min(min(condition_max),0)).^2;

end

% asset grid function

function asset_grid = asset_grid(borrowing_limit,parameters)
    % Asset grid
    I = 100;                            % # of points in the asset grid
    a_grid = power_grid(borrowing_limit,parameters.amax,I);
    n_add = parameters.n_add;
    bl_after = parameters.amin;
    a_append = linspace(bl_after,borrowing_limit,n_add)';
    asset_grid = [a_append(1:n_add-1); a_grid];
end