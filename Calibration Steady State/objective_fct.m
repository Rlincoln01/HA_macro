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
    % retrieve equilibrium objects
    K = ss.capital;
    Y = ss.output;
    T = ss.lump_sum;
    C = ss.consumption;
    w_eq = ss.wage;
    g_final = ss.joint_dist;
    a_marg_dist = ss.wealth_dist;
    r_eq = ss.interest;

    % Target wealth-income ratio
    wealth_income_ratio = K/Y;

    % Target Cash transfers as share of gdp
    transfers_to_gdp = T/Y;

    % obtain income dist
    % productivity grid
    [~,z_grid] = inc_matrix(parameters);

    % other moments of the wealth distribution
    [~,debt_to_gdp,~,wealth_shares,~,~] = wealth_stats(a_marg_dist, ...
        a_grid,z_grid, ...
        Y,"Quintiles",...
        g_final,w_eq,T);
    
    % store simulated moments
    sim_moments.wealth_income_ratio = wealth_income_ratio;
    sim_moments.debt_to_gdp = 100*debt_to_gdp;
    sim_moments.transfers_to_gdp = 100*transfers_to_gdp;
    sim_moments.p0p20 = wealth_shares(1);
    sim_moments.p20p40 = wealth_shares(2);
    sim_moments.p40p60 = wealth_shares(3);
    sim_moments.p60p80 = wealth_shares(4);
    sim_moments.p80p100 = wealth_shares(5);

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
    a_grid = power_grid(borrowing_limit,parameters.amax,parameters.n_a);
    a_append = linspace(parameters.amin,borrowing_limit,parameters.n_add)';
    asset_grid = [a_append(1:parameters.n_add-1); a_grid];
end