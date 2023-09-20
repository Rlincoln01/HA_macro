% ========================================================================
%                   Estimation of Income Process
% ========================================================================
% Description: 
%
% 
%
%
% Author: Rafael Lincoln
% ========================================================================

% =============== Block 1 - Settings and Parametrization =============== %

% Settings of the simulation
settings.n_sim = 5000;              % Number of individual simulations
settings.t_sim = 20;                % Number of years in the simulation
settings.algorithm = 1;             % Algorithm chosen. if 1, global. If 2, local search
settings.multishare = 0.05;         % use x % of the guesses
settings.n_multi = 10000;             % total number of guesses

% Settings of the grid
settings.n_t = 3;                   % # of points in the transitory inc grid
settings.n_p = 11;                  % # of points in the persistent inc grid
settings.pgrid_width = 5.69800327154487/2;       %
settings.tgrid_width = 3.78503497352302/2;        % Transitory grid width
settings.tgrid_par = 0.792311684654520;
settings.pgrid_par = 0.641191531992662;

% settings of the Parameter space
% ORDER: [beta_t,beta_p,sigma_t,sigma_p,lambda_t,lambda_p];
settings.guess_min = [1e-4,1e-4,0.2,0.2,1e-4,1e-4];
settings.guess_max = [1,1,2,2,0.2,0.2];

% Draws of the simulation
seed_t = 1; % transitory shocks
seed_p = 2; % persistent

% Simulate transitory and persistent income shocks
sim_qtrs = floor(settings.t_sim/0.25) + 1;        % Time under simulation (quarters)
rng(seed_t);
draws.t_shocks = rand([settings.n_sim,sim_qtrs]); % transitory

rng(seed_p);
draws.p_shocks = rand([settings.n_sim,sim_qtrs]); % persistent

% Sample moments
sample_moments_file;

% =============== Block 2 - Estimation of the income process =============== %

[par_sol, nfeval] = estimation(settings,draws,sample_moments.ARG,"argentina"); 

% =============== Block 3 - Settings and Parametrization =============== %

% Ex

sigma_t = 1.6249; 
sigma_p = 0.905;
% drift parameters
beta_t = 0.5504;
beta_p = 0.0009;
% poisson arrival rates
lambda_t = 0.0919;
lambda_p = 0.0129;

parameters = [beta_t,beta_p,sigma_t,sigma_p,lambda_t,lambda_p];


yannsim = simulateIncProcess(parameters,settings,draws);

sim_moments = moments(yannsim,settings);

yann = exp(yannsim(:, 1));
d1lyann = yannsim(:, 2) - yannsim(:, 1);
d5lyann = yannsim(:, 5) - yannsim(:, 1);

% Compute relative mean earnings and log of relative mean earnings
yann_relmean = yann / mean(yann);
lyann_relmean = log(yann_relmean);

% Compare with normal dist
y = -5:0.1:5;
mu = 0;
[g_dy1,xi_1] = ksdensity(d1lyann,unique(d1lyann),"Bandwidth",0.1);
[g_dy5,xi_5] = ksdensity(d5lyann,unique(d5lyann),"Bandwidth",0.2);
sigmady1 = sim_moments.std_1yr_log_change;
sigmady5 = sim_moments.std_5yr_log_change;
f_dy1 = normpdf(y,mu,sigmady1);
f_dy5 = normpdf(y,mu,sigmady5);


nbin = 35;

f = figure;
f.Position = [93,44,1202,753];
subplot(2,2,1)
histogram(yann_relmean, nbin, 'Normalization', 'probability', 'FaceColor', 'b', 'EdgeColor', 'k');
tlt = title("$y_{t}$","FontSize",20);
tlt.Interpreter = "latex";
ylabel('Fraction');
ylabel('Fraction');

subplot(2,2,2)
histogram(lyann_relmean, nbin, 'Normalization', 'probability', 'FaceColor', 'b', 'EdgeColor', 'k');
tlt = title("$\ln y_{t}$","FontSize",20);
tlt.Interpreter = "latex";
ylabel('Fraction');

subplot(2,2,3)
% histogram(d1lyann, nbin, 'Normalization', 'pdf', 'FaceColor', 'b', 'EdgeColor', 'k');
p1 = plot(xi_1,g_dy1,"LineWidth",2);
p1.Color = '#000080';
hold on
p2 = plot(y,f_dy1,"LineWidth",1.5,"LineStyle","--");
p2.Color = "#FF3300";
xlim([-2.5,2.5]);
grid on
% histogram(d1lyann, nbin, 'Normalization', 'probability', 'FaceColor', 'r', 'EdgeColor', 'k');
% histfit(d1lyann,nbin,"kernel")
ylabel('Density','interpreter','latex',"FontSize",16);
lgd = legend(["Income Process","$\mathcal{N}(0,\sigma^{2})$"],...
    AutoUpdate="off",Location="northeast");
lgd.Interpreter = "latex";
lgd.FontSize = 11;
tlt = title("$\ln y_{t+1} - \ln y_{t}$","FontSize",20);
tlt.Interpreter = "latex";
st_dev_text = sprintf('St. Dev = %0.2f', sim_moments.std_1yr_log_change);
kurtosis_text = sprintf('Kurtosis = %0.2f', sim_moments.kurt_1yr_log_change);
P9050_text = sprintf('P9050 = %0.2f', sim_moments.P9050_1yr_change);
P5010_text = sprintf('P5010 = %0.2f', sim_moments.P5010_1yr_change);
dim = [0.15 0.15 0.2 0.2];
annotation('textbox',dim,'String',{st_dev_text, kurtosis_text,P9050_text,P5010_text},'FitBoxToText','on', ...
    "Interpreter","latex","FontSize",14,'BackgroundColor','white');


    
% 5 year log change
subplot(2,2,4)

% histogram(d5lyann, nbin, 'Normalization', 'pdf', 'FaceColor', 'b', 'EdgeColor', 'k');
p1 = plot(xi_5,g_dy5,"LineWidth",2);
p1.Color = '#000080';
hold on
p2 = plot(y,f_dy5,"LineWidth",1.5,"LineStyle","--");
p2.Color = "#FF3300";
% histfit(d5lyann,nbin,"kernel")
xlim([-2.5,2.5]);
grid on
ylabel('Density','interpreter','latex',"FontSize",16);
lgd = legend(["Income Process","$\mathcal{N}(0,\sigma^{2})$"],...
    AutoUpdate="off",Location="northeast");
lgd.Interpreter = "latex";
lgd.FontSize = 11;
tlt = title("$\ln y_{t+5} - \ln y_{t}$","FontSize",20);
tlt.Interpreter = "latex";
st_dev_text = sprintf('St. Dev = %0.2f', sim_moments.std_5yr_log_change);
kurtosis_text = sprintf('Kurtosis = %0.2f', sim_moments.kurt_5yr_log_change);
P9050_text = sprintf('P9050 = %0.2f', sim_moments.P9050_5yr_change);
P5010_text = sprintf('P5010 = %0.2f', sim_moments.P5010_5yr_change);
dim = [0.6 0.15 0.2 0.2];
annotation('textbox',dim,'String',{st_dev_text, kurtosis_text,P9050_text,P5010_text},'FitBoxToText','on', ...
    "Interpreter","latex","FontSize",14,'BackgroundColor','white');



