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
specification.disc_rate_heterogeneity = 1;     % Equilibrium with discount rate heterogeneity

% ======================= Parametrization =============================== %

% Income process calibration for US
parameters.sigma_t = 1.74;                     % Variance of transitory shocks
parameters.sigma_p = 1.53;                     % Variance of persistent shocks
parameters.beta_t = 0.761;                     % Persistence of transitory shocks
parameters.beta_p = 0.009;                     % Persistence of permanent shocks
parameters.lambda_p = 0.007;                   % Arrival of transitory shocks
parameters.lambda_t = 0.08;                    % Arrival of permanent shocks


% Income process calibration for Brazil
parameters.sigma_t = 1.6927;                   % Variance of transitory shocks
parameters.sigma_p = 1.1716;                   % Variance of persistent shocks
parameters.beta_t = 0.7910;                    % Persistence of transitory shocks
parameters.beta_p = 0.0019;                    % Persistence of permanent shocks
parameters.lambda_p = 0.0082;                  % Arrival of transitory shocks
parameters.lambda_t = 0.1776;                  % Arrival of permanent shocks


% other parameters
parameters.sigma_omega = 0.5;                  % Variance of the individual fixed effect
parameters.gamma = 1;                          % Relative risk aversion
parameters.alpha = 0.33;                       % Share of capital in production
parameters.dep = 0.025;                        % Depreciation
parameters.r_open_econ = 0.03;                 % Open-economy interest rate
parameters.amax = 600;                         % Maximum asset value
parameters.pgrid_width = 5.69800327154487/2;   % Persistent grid width
parameters.tgrid_width = 3.78503497352302/2;   % Transitory grid width
parameters.tgrid_par = 0.792311684654520;      % Level of non-linearity of the grid
parameters.pgrid_par = 0.641191531992662;      % Level of non-linearity of the grid
parameters.n_t = 3;                            % # of points in the transitory inc grid
parameters.n_p = 11;                           % # of points in the persistent inc grid
parameters.n_rho = 5;                          % # of points on the discount rate grid
parameters.n_fe = 3;                           % # of points in the persistent inc grid
parameters.N_t = 300;                          % # of points in the time grid
parameters.T = 75;                             % termination period of the time grid    


% Calibration
parameters.delta_rho_1 = 0.0033;               % distribution parameters of disc. rate het.
parameters.delta_rho_2 = 0.0020;               % distribution parameters of disc. rate het.
borrowing_limit = -2.5414;                     % Borrowing Limit
parameters.rho = 0.0743;                       % Average discount rate of the economy
parameters.tau = 0.2170;                       % Tax rate

% ============================ Asset grid =============================== %

% Asset grid
I = 100;                            % # of points in the asset grid
a_grid_bef = power_grid(borrowing_limit,parameters.amax,I);
n_add = 21;
bl_after = -5;
a_append = linspace(bl_after,borrowing_limit,n_add)';
a_grid = [a_append(1:n_add-1); a_grid_bef];

% ======================= Call Steady State ============================= %

% Equilibrium distribution, output and prices
tic;
initial_SS = stationary_equilibrium(parameters,specification, ...
                                                   a_grid,borrowing_limit);
toc;

wealth_income_ratio = initial_SS.capital/initial_SS.output;
cash_transfers_to_gdp = 100*initial_SS.lump_sum/initial_SS.output; 

% productivity grid
[lambda_z,z_grid] = inc_matrix(parameters);
n_z = parameters.n_p*parameters.n_t;
dz_tilde = grid_measure(z_grid,n_z); % Define measure for the productivity grid for integration

% Wealth statistics
[gini,debt_to_gdp,htm,wealth_shares,indebted_hhs,lorenz_df] = wealth_stats(initial_SS.wealth_dist, ...
                                                              a_grid,z_grid, ...
                                                              initial_SS.output,"Quintiles",...
                                                              initial_SS.joint_dist, ...
                                                              initial_SS.wage,initial_SS.lump_sum);

% Quarterly MPC
qtr = 1;
windfall_gain = 0.03; % 500 $ in assets
MPCs = Compute_MPCs(qtr,windfall_gain,parameters,specification,0,a_grid,...
                        borrowing_limit, ...
                        initial_SS.interest, ...
                        initial_SS.wage, ...
                        initial_SS.lump_sum);


% Compute average MPC
da_tilde = grid_measure(a_grid,length(a_grid));
avg_mpc = 100*(MPCs.*initial_SS.joint_dist*dz_tilde)'*da_tilde; % aggregate consumption 

% ================= Credit constraint MIT shock ========================= %

% Time grid
N_t = 300; % Final period - by this time, the system should've converged to the new SS
T = 75;    % Number of total quarters ~ 19 yrs
dt = T/N_t; % each dt period corresponds to a quarter

% return to the previous SS
epsilon = 1;
bl_after = borrowing_limit - epsilon;

% Compute transitory shock to the borrowing limit
nu = 0.2;
bl_path = zeros(N_t,1);
% bl_path(1) = borrowing_limit; bl_path(N_t) = bl_after;
bl_path(1) = borrowing_limit; bl_path(N_t) = borrowing_limit;
bl_path(2) = borrowing_limit - epsilon;
bl_after = borrowing_limit;
for n=2:N_t-1
    bl_path(n+1) = bl_path(n) + dt*nu*(bl_after-bl_path(n));
end


% Bang adjustment shock
% a_append = linspace(borrowing_limit,a_min_final,21);
% 
% bl_path = [borrowing_limit, a_min_final*ones(1,N_t-1)];


% Non-linear adjustment path
% nu = 0.2;
% bl_path = zeros(N_t,1);
% bl_path(1) = borrowing_limit; bl_path(N_t) = bl_after;
% for n=1:N_t-1
%     bl_path(n+1) = bl_path(n) + dt*nu*(a_min_final-bl_path(n));
% end
tic;
transition_allocation.levels = transition_path_new(...
    parameters,specification,a_grid,bl_path);
toc;

 
% gini, debt-to-dgp and indebted HHs stats
gini_t = zeros(N_t,1);
debt_to_gdp_t = zeros(N_t,1);
indebted_hhs_t = zeros(N_t,1);

for n=1:N_t
    [gini_t(n),debt_to_gdp_t(n),~,~,indebted_hhs(n),~] = wealth_stats( ...
        transition_allocation.levels.wealth_dist(:,n), ...
        a_grid,z_grid, ...
        transition_allocation.levels.output(n),"Quintiles",...
        transition_allocation.levels.joint_dist(:,:,n), ...
        transition_allocation.levels.wage(n), ...
        transition_allocation.levels.lump_sum(n));
end    

wealth_output_t = (transition_allocation.levels.capital)./transition_allocation.levels.output;

transition_allocation.levels.gini = gini_t;
transition_allocation.levels.debt_to_gdp = 100.*debt_to_gdp_t;
transition_allocation.levels.indebted_hhs = 100.*indebted_hhs_t;
transition_allocation.levels.wealth_output = wealth_output_t;

% Compute % deviations from steady state
% transition_allocation.dev_ss.output = 100.*log(transition_allocation.levels.output./initial_SS.output);
% transition_allocation.dev_ss.consumption = 100.*log(transition_allocation.levels.consumption./initial_SS.consumption);
% transition_allocation.dev_ss.investment = 100.*log(transition_allocation.levels.investment./initial_SS.investment);
% transition_allocation.dev_ss.capital = 100.*log(transition_allocation.levels.capital./initial_SS.capital);
% transition_allocation.dev_ss.interest = 10000.*(transition_allocation.levels.interest - initial_SS.interest); % Basis Points
% transition_allocation.dev_ss.wage = 100.*log(transition_allocation.levels.wage./initial_SS.wage);
% transition_allocation.dev_ss.lump_sum = 100.*log(transition_allocation.levels.consumption./initial_SS.consumption);
% transition_allocation.dev_ss.bl_path = 100.*log(bl_path./bl_path(1));

% Compute % deviations from steady state
transition_allocation.dev_ss.output = 100.*log(transition_allocation.levels.output./transition_allocation.levels.output(1));
transition_allocation.dev_ss.consumption = 100.*log(transition_allocation.levels.consumption./transition_allocation.levels.consumption(1));
transition_allocation.dev_ss.investment = 100.*log(transition_allocation.levels.investment./transition_allocation.levels.investment(1));
transition_allocation.dev_ss.capital = 100.*log(transition_allocation.levels.capital./transition_allocation.levels.capital(1));
transition_allocation.dev_ss.interest = 10000.*(transition_allocation.levels.interest - transition_allocation.levels.interest(1)); % Basis Points
transition_allocation.dev_ss.wage = 100.*log(transition_allocation.levels.wage./transition_allocation.levels.wage(1));
transition_allocation.dev_ss.lump_sum = 100.*log(transition_allocation.levels.lump_sum./transition_allocation.levels.lump_sum(1));
transition_allocation.dev_ss.bl_path = 100.*log(bl_path./bl_path(1));



% Plot the solution
f1= figure(1);
f1.Position = [209,115,929,617];
set(gca,'FontSize',16);
% wealth shares
subplot(3,3,1)
plot(dt*(1:N_t),transition_allocation.dev_ss.output,'LineWidth',2)
grid on
xlim([0,75])        
xlabel('quarters','Interpreter','latex','FontSize',14)
ylabel('\% deviation from S.S.','Interpreter','latex','FontSize',14)
tlt = title("Output - $Y_t$");
tlt.Interpreter = "latex";
tlt.FontSize = 20;
% asset prices
subplot(3,3,2)
plot(dt*(1:N_t),transition_allocation.dev_ss.consumption,'LineWidth',2)
grid on
xlim([0,75])
xlabel('quarters','Interpreter','latex','FontSize',14)
% ylabel('% deviation from S.S.','Interpreter','latex','FontSize',14)
tlt = title("Aggregate Consumption - $C_t$");
tlt.Interpreter = "latex";
tlt.FontSize = 20;
% leverage
subplot(3,3,3)
plot(dt*(1:N_t),transition_allocation.dev_ss.investment,'LineWidth',2)
grid on
xlim([0,75])
xlabel('quarters','Interpreter','latex','FontSize',14)
% ylabel('% deviation from S.S.','Interpreter','latex','FontSize',14)
tlt = title("Aggregate Investment - $I_t$");
tlt.Interpreter = "latex";
tlt.FontSize = 20;
% Volatility of asset prices
subplot(3,3,4)
plot(dt*(1:N_t),transition_allocation.dev_ss.capital,'LineWidth',2)
grid on
xlim([0,75])
xlabel('quarters','Interpreter','latex','FontSize',14)
ylabel('\% deviation from S.S.','Interpreter','latex','FontSize',14)
tlt = title("Aggregate Capital - $K_t$");
tlt.Interpreter = "latex";
tlt.FontSize = 20;
% capital fraction of experts
subplot(3,3,5)
plot(dt*(1:N_t),transition_allocation.dev_ss.interest,'LineWidth',2)
grid on
xlim([0,75])
xlabel('quarters','Interpreter','latex','FontSize',14)
% ylabel('% deviation from S.S.','Interpreter','latex','FontSize',14)
tlt = title("Real interest rate - $r_t$ (b.p.)");
tlt.Interpreter = "latex";
tlt.FontSize = 20;
% risk borne by experts
subplot(3,3,6)
plot(dt*(1:N_t),transition_allocation.dev_ss.wage,'LineWidth',2)
grid on
xlim([0,75])
xlabel('quarters','Interpreter','latex','FontSize',14)
% ylabel('% deviation from S.S.','Interpreter','latex','FontSize',14)
tlt = title("Wages - $w_t$");
tlt.Interpreter = "latex";
tlt.FontSize = 20;
% Borrowing Limit path
subplot(3,3,7)
plot(dt*(1:N_t),transition_allocation.dev_ss.bl_path,'LineWidth',2)
grid on
xlim([0,75])
xlabel('quarters','Interpreter','latex','FontSize',14)
ylabel('\% deviation from S.S.','Interpreter','latex','FontSize',14)
tlt = title("Borrowing Limit $\underbar{a}$");
tlt.Interpreter = "latex";
tlt.FontSize = 20;
% Gini
subplot(3,3,8)
plot(dt*(1:N_t),transition_allocation.levels.gini,'LineWidth',2)
grid on
xlim([0,75])
xlabel('quarters','Interpreter','latex','FontSize',14)
% ylabel('% deviation from S.S.','Interpreter','latex','FontSize',14)
tlt = title("Wealth Gini");
tlt.Interpreter = "latex";
tlt.FontSize = 20;
% Debt-to-GDP
subplot(3,3,9)
plot(dt*(1:N_t),transition_allocation.levels.debt_to_gdp,'LineWidth',2)
grid on
xlim([0,75])
xlabel('quarters','Interpreter','latex','FontSize',14)
% ylabel('% deviation from S.S.','Interpreter','latex','FontSize',14)
tlt = title("Debt to GDP - $\frac{D_t}{Y_t}$");
tlt.Interpreter = "latex";
tlt.FontSize = 20;
sgtitle('Credit Shock','Interpreter','latex','FontSize',20);

% ================= IRF MIT-shock borrowing constraint ========================= %

% general equilibrium consumption
C_t_GE = transition_allocation.levels.consumption;

% Partial Eq. paths for theta = \underline{a}
Gamma_t.PE_credit.wages_t = initial_SS.wage.*ones(1,N_t);
Gamma_t.PE_credit.bl_t = bl_path;
Gamma_t.PE_credit.interest_t = initial_SS.interest.*ones(1,N_t);
Gamma_t.PE_credit.lump_sum_t = initial_SS.lump_sum.*ones(1,N_t);

C_t_PE_credit = PE_decomposition(parameters,specification,Gamma_t.PE_credit,a_grid,initial_SS);

% Partial Eq. paths for theta = w_t
Gamma_t.PE_wages.wages_t = transition_allocation.levels.wage;
Gamma_t.PE_wages.bl_t = bl_path(1).*ones(1,N_t);
Gamma_t.PE_wages.interest_t = initial_SS.interest.*ones(1,N_t);
Gamma_t.PE_wages.lump_sum_t = initial_SS.lump_sum.*ones(1,N_t);

C_t_PE_wages = PE_decomposition(parameters,specification,Gamma_t.PE_wages,a_grid,initial_SS);

% Partial Eq. paths for theta = r_t
Gamma_t.PE_interest.wages_t = initial_SS.wage.*ones(1,N_t);
Gamma_t.PE_interest.bl_t = bl_path(1).*ones(1,N_t);
Gamma_t.PE_interest.interest_t = transition_allocation.levels.interest;
Gamma_t.PE_interest.lump_sum_t = initial_SS.lump_sum.*ones(1,N_t);

C_t_PE_interest = PE_decomposition(parameters,specification,Gamma_t.PE_interest,a_grid,initial_SS);

% Partial Eq. paths for theta = lump_sum_t
Gamma_t.PE_lump_sum.wages_t = initial_SS.wage.*ones(1,N_t);
Gamma_t.PE_lump_sum.bl_t = bl_path(1).*ones(1,N_t);
Gamma_t.PE_lump_sum.interest_t = initial_SS.interest.*ones(1,N_t);
Gamma_t.PE_lump_sum.lump_sum_t = transition_allocation.levels.lump_sum;

C_t_PE_lump_sum = PE_decomposition(parameters,specification,Gamma_t.PE_lump_sum,a_grid,initial_SS);

% ==== Compute Deviations from SS ==== %

dC_t_GE = 100.*log(C_t_GE./initial_SS.consumption);
dC_t_PE_lump_sum = 100.*log(C_t_PE_lump_sum./initial_SS.consumption);
dC_t_PE_interest = 100.*log(C_t_PE_interest./initial_SS.consumption);
dC_t_PE_wages = 100.*log(C_t_PE_wages./initial_SS.consumption);
dC_t_PE_credit = 100.*log(C_t_PE_credit./initial_SS.consumption);

f = figure(2);
plot(dt*(1:(N_t-1)),dC_t_GE(2:end),LineWidth=2,LineStyle="-")
hold on
plot(dt*(1:(N_t-1)),dC_t_PE_credit(2:end),LineWidth=2,LineStyle="-.")
hold on
plot(dt*(1:(N_t-1)),dC_t_PE_wages(2:end) + dC_t_PE_lump_sum(2:end),LineWidth=2,LineStyle=":")
hold on
plot(dt*(1:(N_t-1)),dC_t_PE_interest(2:end),LineWidth=2,LineStyle="--")
grid on
yline(0)
legend("Total Effect","Direct: $\underline{a}_{t}$", ...
    "Indirect: $w_t + T_t$","Indirect: $r_t$","interpreter","latex", ...
    "fontsize",20)
xlabel('Quarters','Interpreter','latex','FontSize',18)
ylabel('\% deviation from S.S.','Interpreter','latex','FontSize',18)
tlt = title({"Consumption Decomposition: Partial ($\underline{a}$)", ...
    "and General Equilibrium ($\left\{ r_t,w_t,T_t\right\}$) effects"});
tlt.Interpreter = "latex";
tlt.FontSize = 20;
saveas(f,fullfile(filepath,"con_decomposition"),"epsc");


% ================= Computing elasticities ========================= %

% compute aggregate variables' deviation from SS
dtot.output = transition_allocation.levels.output - initial_SS.output;
dtot.consumption = transition_allocation.levels.consumption - initial_SS.consumption;
dtot.investment = transition_allocation.levels.investment - initial_SS.investment;
dtot.C_PE_credit = C_t_PE_credit - initial_SS.consumption;
dtot.C_PE_lump_sum = C_t_PE_lump_sum - initial_SS.consumption;
dtot.C_PE_interest = C_t_PE_interest - initial_SS.consumption;
dtot.C_PE_wages = C_t_PE_wages - initial_SS.consumption;

% Generate tables of n quarters of elasticity of variables

% periods: (dt = 1 qtr)
tset{1} = [1:3]+1;   % 3 quarters
tset{2} = [1:6]+1;   % 6 quarters
tset{3} = [1:9]+1;   % 9 quarters
tset{4} = [1:12]+1;  % 12 quarters
tset{5} = [1:75] + 1;

% % periods for \underline{a}
% for i = 1:4
%     tsetA{i} = tset{i}+1;
% end

tab = zeros(13,5);
elastdenom = zeros(1,5);

for col = 1:5
    t_set = tset{col};
    % t_setA = tsetA{col};

    % borrowing limit change
    row = 1;
    tab(row,col) = sum((bl_path(t_set)-bl_path(1))./bl_path(1))/length(t_set);
    elastdenom(1,col) = tab(row,col);
    % Real interest rate change
    row = row + 1;
    tab(row,col) = sum(transition_allocation.interest(t_set))/length(t_set) - initial_SS.interest;
    % wage changes
    row = row + 1;
    tab(row,col) = sum(transition_allocation.wage(t_set))/length(t_set) - initial_SS.wage;
    % output change and elasticity
    row = row + 1;
    tab(row,col) = sum(dtot.output(t_set)./initial_SS.output)/length(t_set);
    row = row + 1;
    tab(row,col) = tab(row-1,col) ./ elastdenom(1,col) ;
    % consumption change and elasticity
    row = row + 1;
    tab(row,col) = sum(dtot.consumption(t_set)./initial_SS.consumption)/length(t_set);
    row = row + 1;
    tab(row,col) = tab(row-1,col) ./ elastdenom(1,col) ;
    % investment change and elasticity
    row = row + 1;
    tab(row,col) = sum(dtot.investment(t_set)./initial_SS.investment)/length(t_set);
    row = row + 1;
    tab(row,col) = tab(row-1,col) ./ elastdenom(1,col) ;

    % decomposition of the elasticity of aggregate consumption by direct
    % and indirect effs

    % direct effects - credit
    row = row+1;
    tab(row,col) = sum(dtot.C_PE_credit(t_set))./sum(dtot.consumption(t_set));
    % indirect effects - interest
    row = row+1;
    tab(row,col) = sum(dtot.C_PE_interest(t_set))./sum(dtot.consumption(t_set));
    % indirect effects - wages
    row = row + 1;
    tab(row,col) = sum(dtot.C_PE_wages(t_set))./sum(dtot.consumption(t_set));
    % indirect effects - lump_sum
    row = row + 1;
    tab(row,col) = sum(dtot.C_PE_lump_sum(t_set))./sum(dtot.consumption(t_set));
    
end    

format long;
disp(tab);
 

% ================= Lorenz Curve ========================= %

% Lorenz curve for the SCF data 2019 - NyFED
scf_lorenz = 100.*[0.0002 -0.0002; ...
                0.1027 	-0.0054;...
                0.2021 	-0.0051;...
                0.3051 	-0.0031;...
                0.4203 	0.0048;...
                0.5429 	0.0230;...
                0.6640 	0.0565;...
                0.7757 	0.1098;...
                0.8682 	0.1925;...
                0.9449 	0.3351;...
                0.9861 	0.5792;...
                0.9993 	0.8810;...
                1 	1];

% MPC and Wealth distributions
asset_to_usd = 16.727; % Mapping Thousands of USD to a=1;
mean_wealth = initial_SS.capital;
wage = initial_SS.wage;
wealth_dist_cut = initial_SS.wealth_dist(n_add:end);
asset_usd_grid = asset_to_usd.*a_grid_bef;
wealth_prob_dist = (initial_SS.wealth_dist).*da_tilde;

% Integrate the income part to obtain the marginal dist of MPCs
MPC_wealth_dist = 100.*(MPCs.*initial_SS.joint_dist*dz_tilde);
MPC_wealth_dist = MPC_wealth_dist(n_add:end);

plot(scf_lorenz(:,1),scf_lorenz(:,2))
hold on
plot(lorenz_df(:,1),lorenz_df(:,2))

% ======================== Lorenz Curve ======================== %
f = figure(4);
f.Position = [93,336,1061,461];
tiledlayout(1,2, 'TileSpacing', 'compact'); 
nexttile 
linha = 0:100;
% f= figure(3);
% f.Position = [151,208,700,585];
set(gca,'FontSize',20)
p1 = plot(scf_lorenz(:,1),scf_lorenz(:,2),"LineWidth",2);
p1.Color = '#000080';
hold on
p2 = plot(lorenz_df(:,1),lorenz_df(:,2),"LineWidth",2);
p2.Color = '#4682B4';
hold on
p4 = plot(linha,linha,"--","LineWidth",1);
p4.Color = "black";
grid on;
yline(0);
% Ticks
xticks(0:10:100);
yticks(-10:10:100);
ylim([-5,100]);
xlim([0,100]);
% Text
text(50,55,"Line of Equality - $45^{\circ}$","rotation",40,"Interpreter","Latex", ...
    "HorizontalAlignment","center","FontSize",14)
% Labels
xlabel("Cumulative \% of Population by Percentile",'interpreter','latex',"FontSize",16)
ylabel('Cumulative \% of Net Wealth','interpreter','latex', ...
    "FontSize",16)
% Title
tlt = title("Net Wealth - Lorenz Curves");
tlt.Interpreter = "latex";
tlt.FontSize = 20;
% Legend
lgd = legend(["SCF (US) Data - 2019", ...
                "Model"], ...
                AutoUpdate="off", ...
                Location="northwest");
lgd.Interpreter = "latex";
lgd.FontSize = 18;
% ================= Wealth & MPC Distribution ======================== %
nexttile
% f = figure;
% f.Position = [93,234,750,563];
set(gca,'FontSize',20)
yyaxis left
bar(asset_usd_grid,wealth_dist_cut,1,"histc");
ylabel("Density",'Interpreter','latex','FontSize',20)
sh  = findall(gcf,'marker','*'); delete(sh);
yyaxis right
plot(asset_usd_grid,MPC_wealth_dist,linewidth = 3, ...
    LineStyle="-.");
xlim([-50,200])
grid on
xlabel('US\$ Thousands (2022)','Interpreter','latex','FontSize',20)
ylab = ylabel('\% Consumed out of 500\$ in a quarter','FontSize',20,'interpreter','latex', ...
    'Rotation',90);
x1= xline(asset_usd_grid(1),"k--", ...
    {"$\underline{a} = $" + sprintf(" %.3f",asset_usd_grid(1)) + "\$", ...
    "$P(a = \underline{a}) = $ " + sprintf(" %.2f",max(wealth_prob_dist))});
x1.LineWidth = 1;
x1.Interpreter = "latex";
x1.FontSize = 18;
x1.LabelOrientation = "horizontal";
tlt = title("Quarterly MPCs out of 500\$", "FontSize",26);
tlt.Interpreter = "latex";
saveas(f,fullfile(filepath,"wealth_distribution_data_model"),"epsc");


