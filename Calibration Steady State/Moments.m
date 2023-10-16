% ========================================================================
%                   Steady State Moments Function
% ========================================================================
% Description: Generates the steady-state moments targeted in the data
%
% 
%
%
% Author: Rafael Lincoln
% ========================================================================

function sim_moments = Moments(ss,parameters,a_grid)

% retrieve equilibrium objects
K = ss.capital;
Y = ss.output;
T = ss.lump_sum;
w_eq = ss.wage;
g_final = ss.joint_dist;
a_marg_dist = ss.wealth_dist;

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
sim_moments.p0p20 = 100*wealth_shares(1);
sim_moments.p20p40 = 100*wealth_shares(2);
sim_moments.p40p60 = 100*wealth_shares(3);
sim_moments.p60p80 = 100*wealth_shares(4);
sim_moments.p80p100 = 100*wealth_shares(5);

end