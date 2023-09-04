% =========================================================================
% Main File:            
%
%
% • Parametrization:
% • Model specification:
%
%
%
% Author: Rafael Lincoln
% =========================================================================


% ======================= Model Specification =========================== %

specification.ind_fixed_effects = 0;           % Compute individual fixed effects
specification.closed_economy = 1;              % Compute closed economy interest rate
specification.disc_rate_heterogeneity = 0;     % Equilibrium with discount rate heterogeneity

% ======================= Parametrization =============================== %

parameters.sigma_t = 1.74;                     % Variance of transitory shocks
parameters.sigma_p = 1.53;                     % Variance of persistent shocks
parameters.beta_t = 0.761;                     % Persistence of transitory shocks
parameters.beta_p = 0.009;                     % Persistence of permanent shocks
parameters.lambda_p = 0.007;                   % Arrival of transitory shocks
parameters.lambda_t = 0.08;                    % Arrival of permanent shocks
parameters.sigma_omega = 0.5;                  % Variance of the individual fixed effect
parameters.gamma = 2;                          % Relative risk aversion
parameters.rho = 0.05;                         % discount rate of households
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


% ============================ Asset grid =============================== %

% Asset grid
borrowing_limit = -1;               % Borrowing limit for liquid wealth
I = 100;                            % # of points in the asset grid
a_grid = power_grid(borrowing_limit,parameters.amax,I);
n_add = 21;
bl_after = -5;
a_append = linspace(bl_after,borrowing_limit,n_add)';
a_grid = [a_append(1:n_add-1); a_grid];

% ======================= Call Steady State ============================= %

% Equilibrium distribution, output and prices
[Y,C,S,K_eq,L,lump_sum,r_eq,...
w_eq,g_final,a_marg_dist] = stationary_equilibrium(parameters,specification, ...
                                                   a_grid,borrowing_limit);


plot(a_grid./K_eq,a_marg_dist,"LineWidth",2)
xlim([-1,2])
xlabel('Liquid wealth w.r.t mean wealth, $a$','Interpreter','latex','FontSize',14)
ylabel("Density",'Interpreter','latex','FontSize',14)


% productivity grid
[lambda_z,z_grid] = inc_matrix(parameters);
n_z = parameters.n_p*parameters.n_t;
dz_tilde = grid_measure(z_grid,n_z); % Define measure for the productivity grid for integration

% Wealth statistics
[gini,debt_to_gdp,htm,wealth_shares,indebted_hhs,lorenz_df] = wealth_stats(a_marg_dist, ...
                                                              a_grid,z_grid, ...
                                                              Y,"Quintiles",...
                                                              g_final,w_eq,lump_sum);

% Quarterly MPC
qtr = 1;
windfall_gain = 1/50;
MPCs = Compute_MPCs(qtr,windfall_gain,parameters,specification,0,a_grid,...
                        borrowing_limit,r_eq,w_eq,lump_sum);


% Compute average MPC
da_tilde = grid_measure(a_grid,length(a_grid));
avg_mpc = 100*(MPCs.*g_final*dz_tilde)'*da_tilde; % aggregate consumption 

% ================= Credit constraint MIT shock ========================= %

% Time grid
N_t = 300; % Final period - by this time, the system should've converged to the new SS
T = 75;
dt = T/N_t;

borrowing_limit = -1;
a_min_final = -2;


% Bang adjustment shock
a_append = linspace(borrowing_limit,a_min_final,21);

bl_path = [borrowing_limit, a_min_final*ones(1,N_t-1)];


% Non-linear adjustment path
nu = 0.2;
bl_path = zeros(N_t,1);
bl_path(1) = borrowing_limit; bl_path(N_t) = bl_after;
for n=1:N_t-1
    bl_path(n+1) = bl_path(n) + dt*nu*(a_min_final-bl_path(n));
end


[Y_t,C_t,I_t,K_t,gg_t,a_marg_t,r_t,w_t,lump_sum_t] = transition_path(...
    parameters,specification,a_grid,bl_path);

% gini, debt-to-dgp and indebted HHs stats
gini_t = zeros(N_t,1);
debt_to_gdp_t = zeros(N_t,1);
indebted_hhs_t = zeros(N_t,1);

for n=1:N_t
    [gini_t(n),debt_to_gdp_t(n),~,~,indebted_hhs(n),~] = wealth_stats(cell2mat(a_marg_t(n)), ...
                                                              a_grid,z_grid, ...
                                                              Y_t(n),"Quintiles",...
                                                              gg_t(:,:,n),w_t(n),lump_sum_t(n));
end    

% [Y_t,C_t,I_t,K_t,gg_t,a_marg_t,r_t,w_t,lump_sum_t] = Copy_of_transition_path(...
%     a_grid,bl_path,permanent_income,transitory_income,closed_economy,0);


wealth_output_t = K_t./Y_t;

% Plot the solution
f1= figure(1);
f1.Position = [209,115,929,617];
set(gca,'FontSize',16);
% wealth shares
subplot(3,3,1)
plot(dt*(1:N_t),Y_t,'LineWidth',2)
xlim([0,75])    
xlabel('t','Interpreter','latex','FontSize',14)
ylabel('$Y_t$','Interpreter','latex','FontSize',14)
tlt = title("Output");
tlt.Interpreter = "latex";
tlt.FontSize = 20;
% asset prices
subplot(3,3,2)
plot(dt*(1:N_t),C_t,'LineWidth',2)
xlim([0,75])
xlabel('t','Interpreter','latex','FontSize',14)
ylabel('$C_t$','Interpreter','latex','FontSize',14)
tlt = title("Aggregate Consumption");
tlt.Interpreter = "latex";
tlt.FontSize = 20;
% leverage
subplot(3,3,3)
plot(dt*(1:N_t),I_t,'LineWidth',2)
xlim([0,75])
xlabel('t','Interpreter','latex','FontSize',14)
ylabel('$I_t$','Interpreter','latex','FontSize',14)
tlt = title("Aggregate Investment");
tlt.Interpreter = "latex";
tlt.FontSize = 20;
% Volatility of asset prices
subplot(3,3,4)
plot(dt*(1:N_t),K_t,'LineWidth',2)
xlim([0,75])
xlabel('t','Interpreter','latex','FontSize',14)
ylabel('$K_t$','Interpreter','latex','FontSize',14)
tlt = title("Aggregate Capital");
tlt.Interpreter = "latex";
tlt.FontSize = 20;
% capital fraction of experts
subplot(3,3,5)
plot(dt*(1:N_t),r_t,'LineWidth',2)
xlim([0,75])
xlabel('t','Interpreter','latex','FontSize',14)
ylabel('$r_t$','Interpreter','latex','FontSize',14)
tlt = title("Real interest rate");
tlt.Interpreter = "latex";
tlt.FontSize = 20;
% risk borne by experts
subplot(3,3,6)
plot(dt*(1:N_t),w_t,'LineWidth',2)
xlim([0,75])
xlabel('t','Interpreter','latex','FontSize',14)
ylabel('$w_t$','Interpreter','latex','FontSize',14)
tlt = title("Wages");
tlt.Interpreter = "latex";
tlt.FontSize = 20;
% Borrowing Limit path
subplot(3,3,7)
plot(dt*(1:N_t),bl_path,'LineWidth',2)
xlim([0,75])
xlabel('t','Interpreter','latex','FontSize',14)
ylabel('$\underbar{a}$','Interpreter','latex','FontSize',14)
tlt = title("Borrowing Limit $\underbar{a}$");
tlt.Interpreter = "latex";
tlt.FontSize = 20;
% Gini
subplot(3,3,8)
plot(dt*(1:N_t),gini_t,'LineWidth',2)
xlim([0,75])
xlabel('t','Interpreter','latex','FontSize',14)
ylabel('$gini_t$','Interpreter','latex','FontSize',14)
tlt = title("Wealth Gini");
tlt.Interpreter = "latex";
tlt.FontSize = 20;
% Debt-to-GDP
subplot(3,3,9)
plot(dt*(1:N_t),debt_to_gdp_t,'LineWidth',2)
xlim([0,75])
xlabel('t','Interpreter','latex','FontSize',14)
ylabel('$\frac{D_t}{Y_t}$','Interpreter','latex','FontSize',14)
tlt = title("Debt to GDP");
tlt.Interpreter = "latex";
tlt.FontSize = 20;
sgtitle('Credit Shock','Interpreter','latex','FontSize',20);







