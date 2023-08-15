% ========================================================================
%                   Objective function
% ========================================================================
% Description: Runs the simulation, computes the simulated moments
%              and calculates the loss function
%
% Inputs:
%   - parameters: set of parameters that govern the income process
%
%
% Author: Rafael Lincoln
% ========================================================================

function loss = objective_fct(parameters,settings,draws,sample_moments)

% Penalty function => fminsearch has no bound restrictions
penalty = penalty_fct(parameters,settings); %function checks whether param bounds have been violated

% ============ Simulate and compute moments ============ %
if penalty == 0
    yannsim = simulateIncProcess(parameters,settings,draws);
    sim_moments = moments(yannsim,settings);
    loss = loss_fct(sample_moments,sim_moments);
else
    loss = 10e6;
end

end


function penalty = penalty_fct(parameters,settings)
    % Preallocation
    penalty= 0;
    % If a parameter is below the minimum values
    condition_min = parameters >= settings.guess_min;
    penalty = penalty + 10000*(min(min(condition_min),0)).^2;
    % If a parameter is above the maximum values
    condition_max = parameters <= settings.guess_max;
    penalty = penalty + 10000*(min(min(condition_max),0)).^2;

end