%=========================================================================
%                 Wealth Statistics
%=========================================================================
% Calculate wealth statistics used in the model
%
% Inputs:
%   - a_grid: Wealth grid used in calculations
%   - z_grid: Productivity grid used in calculations
%   - Y: Output
%   - a_marg_dist: Marginal Wealth distribution
%   - wealth_shares: Two options for aggregation
%       1) "Quintiles" - 0p20,20p40,…,80p100
%       2) "Bottom50-Top1" - bottom 50, 50p90, 90p99, top1
%   - g: Joint productivity-wealth distribution
%   - w: wages
%
%
% Output:  
%   - lorenz_df: table with the x and y values of the lorenz curve
%   - gini: Gini index of net-wealth
%   - wealth_shares: wealth shares according to the group selected
%   - debt_to_gdp: Debt-to-gdp
%   - htm: share of households in the distribution which are HtM
%   - indebted_hhs: Share of households with negative liquid wealth (indebted) 
%
%=========================================================================

function [gini,debt_to_gdp,htm,wealth_shares,indebted_hhs,lorenz_df] = wealth_stats(a_marg_dist,a_grid,z_grid,Y,ws_composition,...
                                        g,w,lump_sum)


% 1) Preliminaries

n_a = length(a_grid);
n_z = length(z_grid);

% asset and productivity grid measure
da_tilde = grid_measure(a_grid,n_a);
dz_tilde = grid_measure(z_grid,n_z); % Define measure for the productivity grid for integration

%=========================================================================
%           Block 1: Lorenz curve and Gini
%=========================================================================

% 1st - Find lorenz curve
% Average Savings/Assets of distribution
integrand = a_grid.*a_marg_dist.*da_tilde;
total = sum(integrand);

% Cumulative % of Wealth
cum_wealth = 100.*(cumsum(integrand)./abs(total));

% Cumulative % of population - wealth
cum_pop_wealth = 100*cumsum(a_marg_dist.*da_tilde);

lorenz_df = [0 0;cum_pop_wealth cum_wealth];

% Distance of the Lorenz curve to the equality line
dist = lorenz_df(:,1) - lorenz_df(:,2);

% 2nd - Calculate gini from lorenz curve
% Area over which I calculate gini,
x_axis = max(lorenz_df(:,1),dist);

% Measure over which I calculate the integral
measure = [0; diff(lorenz_df(:,1))];

% Gini Index
gini = sum(dist.*measure)/sum(x_axis.*measure);


%=========================================================================
%           Block 2: Wealth shares and debt by gdp
%=========================================================================

% - Wealth Shares - %
% CDF
G_a = cumsum(a_marg_dist.*da_tilde);

p = linspace(0,1,100);
for pp=1:100
    [val index(pp)] = min(abs(G_a - p(pp)));
end

if ws_composition == "Quintiles"
    Q1 = sum(integrand(1:index(20)-1))/total;
    Q2 = sum(integrand(index(20):index(40)-1))/total;
    Q3 = sum(integrand(index(40):index(60)-1))/total;
    Q4 = sum(integrand(index(60):index(80)-1))/total;
    Q5 = sum(integrand(index(80):n_a))/total;

    wealth_shares = [Q1,Q2,Q3,Q4,Q5];
elseif ws_composition == "Bottom50-Top1"
    bottom50 = sum(integrand(1:index(50)-1))/total;
    next40 = sum(integrand(index(50):index(90)-1))/total;
    next9 = sum(integrand(index(90):index(99)-1))/total;
    top1 = sum(integrand(index(99):n_a))/total;

    wealth_shares = [bottom50,next40,next9,top1];
else
    error("Wealth shares not specified correctly. Choose between 'Quintiles' and 'Bottom50-Top1'")
end    

% - Debt to gdp - %
debt = sum(min(integrand,0));
debt_to_gdp = -debt/Y;

%=========================================================================
%           Block 3: Share of Hand-to-Mouth households
%=========================================================================

zz = ones(n_a,1)*z_grid';                 % Matrix IxM of productivity
aa = a_grid*ones(1,n_z);                  % Matrix IxM of assets

% obtain the matrix of quarterly income of households
yy = (w.*exp(zz)+ lump_sum);

% which households are HtM?
HtM_map = aa < yy./2;

% distribution of HtM households
g_htm = HtM_map.*g; 

% calculate share of HtM
htm_marginal = g_htm*dz_tilde; % marginal htm dist (wealth)
htm = sum(htm_marginal.*da_tilde);


% %% ========================= Indebted and on the bl HH's ================ %%
% 
% % Households on the borrowing limit
% hh_on_bl = G_a(1);
% 
% Households with a < 0
if min(a_grid) == 0
    indebted_hhs = 0;
else
    i = find(a_grid<0);
    i = i(end);
    indebted_hhs = G_a(i);
end    
% 

end