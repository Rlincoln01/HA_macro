%=========================================================================
% Macro 2 - Incomplete markets and heterogenous agents
% Teacher: Fernando Mendo
% Semester: 2022.2
% Rafael Lincoln - 2212890
%
% Exercise 2 - Continuous time Bewley Model with Exogenous Fiscal Policy
%=========================================================================


%=========================================================================
%                   Parameters and calibration
%=========================================================================

% Poisson Intensities
parameters.lambda_12 = 1.5;
parameters.lambda_13 = 0.5;
parameters.lambda_21 = 0.75;
parameters.lambda_23 = 0.75;
parameters.lambda_31 = 0.5;
parameters.lambda_32 = 0.75;
% Income states
y_1 = 0.1; y_2 = 0.15; y_3 = 0.2;
y = [y_1,y_2,y_3];
parameters.income_states = y;
% Tax rate
tau_0 = 0.2;
parameters.borrow_lim = -0.15;
% Preferences
parameters.gamma = 2;
parameters.rho = .05;
% asset grid parametrization
parameters.n_a = 500;
parameters.amax = 3;

% Asset grid (for plotting reasons)
I=parameters.n_a; % size of asset grid
amin = parameters.borrow_lim;
amax = parameters.amax;
a = linspace(amin,amax,I)'; % asset grid
da = (amax-amin)/(I-1); % step on asset grid

%=========================================================================
%                   Stationary distribution of income process
%=========================================================================

lambda_12 = parameters.lambda_12;
lambda_13 = parameters.lambda_13;
lambda_21 = parameters.lambda_21;
lambda_23 = parameters.lambda_23;
lambda_31 = parameters.lambda_31;
lambda_32 = parameters.lambda_32;

% transition matrix
T = [-(lambda_12 + lambda_13) lambda_12 lambda_13;
    lambda_21 -(lambda_21 + lambda_23) lambda_23;
    lambda_31 lambda_32 -(lambda_31 + lambda_32)];

y = parameters.income_states;

ydist = stationary_dist(T,3);

%=========================================================================
%                   Question c) - Partial Eq.
%=========================================================================

r = 0.035;

[V,con,sav,A] = HJB(parameters,r,tau_0);
[gg,g] = KF(parameters,A);

figure(2)
set(gca,'FontSize',16)
plot(a,g,"LineWidth",3)
xlim([amin-0.01,0.5])
x= xline(parameters.borrow_lim,"k--","$a = \bar{a}$");
grid;
x.LineWidth = 2;
x.Interpreter = "latex";
x.FontSize = 14;
lgd = legend(["$y_{1}$","$y_{2}$","$y_{3}$"],AutoUpdate="off");
lgd.Interpreter = "latex";
xlabel('Wealth, $a$','interpreter','latex')
ylabel('Densities, $g_i(a)$','interpreter','latex')
tlt = title("Equilibrium Densities");
tlt.Interpreter = "latex";
tlt.FontSize = 14;

figure(3)
set(gca,'FontSize',16)
plot(a,con,"LineWidth",2)
xlim([amin,1])
lgd = legend(["$y_{1}$","$y_{2}$","$y_{3}$"],AutoUpdate="off");
lgd.Interpreter = "latex";
xlabel('Wealth, $a$','interpreter','latex')
ylabel('Consumption, $c_i(a)$','interpreter','latex')
tlt = title("Consumption policy functions");
tlt.Interpreter = "latex";
tlt.FontSize = 14;

figure(4)
set(gca,'FontSize',16)
plot(a,sav,"LineWidth",2)
xlim([amin,1])
lgd = legend(["$y_{1}$","$y_{2}$","$y_{3}$"],AutoUpdate="off");
lgd.Interpreter = "latex";
xlabel('Wealth, $a$','interpreter','latex')
ylabel('Savings, $s_i(a)$','interpreter','latex')
tlt = title("Savings policy functions");
tlt.Interpreter = "latex";
tlt.FontSize = 14;

%=========================================================================
%                   Question d) - Asset Supply & Demand
%=========================================================================

% == grid of interest rates == %
rmin = 0.0001;
rmax = parameters.rho-rmin;
nr = 100;
rgrid = linspace(rmin,rmax,nr);


% == Bond Demand == %
S_d = zeros(1,nr);

disp("Calculating net asset supply in the economy, please wait...")
for ir = 1:nr
    [V,con,sav,A] = HJB(parameters,rgrid(ir),tau_0);
    [gg,g] = KF(parameters,A);
    S_d(1,ir) = g(:,1)'*a*da + g(:,2)'*a*da + g(:,3)'*a*da;
end
disp("Done")


% Bond Supply - test1

taxes = tau_0*(y*ydist);

B_d = @(r) taxes./r;

B_dem = B_d(rgrid);

figure(5)
plot(B_dem,rgrid,"r-", ...
    S_d,rgrid,"b-","LineWidth",3)
grid;
xlim([-0.2,1])
ylim([0.04,rho+0.0005]);
yl = yline(rho,"k--","$r = \rho$");
yl.Interpreter = "latex";
yl.LineWidth = 2;
x= xline(borrow_lim,"k--","$a = \bar{a}$");
x.LineWidth = 2;
x.Interpreter = "latex";
x.FontSize = 14;
x.LabelVerticalAlignment = 'middle';
x.LabelHorizontalAlignment = 'center';
lgd = legend(["$B(r)$","$S(r)$"],AutoUpdate="off");
lgd.Location = "southeast";
lgd.Interpreter = "latex";
tlt = title("Bond Supply $S(r)$ and Demand $B(r)$");
tlt.Interpreter = "latex";
tlt.FontSize = 14;
xlbl = xlabel("$B$");
ylbl = ylabel("$r$");
xlbl.Interpreter = "latex";
ylbl.Interpreter = "latex";

%=========================================================================
%                   Question e) - Stationary Equilibrium
%=========================================================================

[B_eq,r_eq] = stationary_equilibrium(parameters,ydist,tau_0);












