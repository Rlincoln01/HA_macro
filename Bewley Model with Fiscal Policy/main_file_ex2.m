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

% == Household == %
% Poisson Intensities
lambda_12 = 1.5; lambda_13 = 0.5; lambda_21 = 0.75;
lambda_23 = 0.75; lambda_31 = 0.5; lambda_32 = 0.75;

% Income states
y_1 = 0.1; y_2 = 0.15; y_3 = 0.2;
y = [y_1,y_2,y_3];

% tax
tau_0 = 0.2;

% borrowing constraint
borrow_lim = -0.15;

% Preferences
gamma = 2; % IES
rho = .05; % beta ~ 0.99
u = @(c) (c.^(1-gamma))./(1-gamma); % CRRA preferences

% == Asset & Income grid == %
I=500; % size of asset grid
amin = borrow_lim;
amax = 3;
a = linspace(amin,amax,I)'; % asset grid
da = (amax-amin)/(I-1); % step on asset grid
y = [y_1,y_2,y_3]; %change
display= 0;

% ???
aa = [a,a,a]; 
yy = ones(I,1)*y;


%=========================================================================
%                   Question a)
%=========================================================================

T = [0 lambda_12/(lambda_12 + lambda_13) lambda_13/(lambda_12 + lambda_13);
    lambda_21/(lambda_21 + lambda_23) 0 lambda_21/(lambda_21 + lambda_23);
    lambda_31/(lambda_31 + lambda_32) lambda_32/(lambda_31 + lambda_32) 0];

T = [(lambda_12+lambda_13) (-lambda_21) (-lambda_31);
    (-lambda_12) (lambda_21+lambda_23) (-lambda_32);
    1 1 1];

ydist = inv(T)*[0 0 1]';


%=========================================================================
%                   Question c) - Partial Eq.
%=========================================================================

r = 0.035;

[V,con,sav,A] = HJB_solve_new(r);
[gg,g] = KF_solve_new(A);

figure(2)
set(gca,'FontSize',16)
plot(a,g,"LineWidth",3)
xlim([amin-0.01,0.5])
x= xline(borrow_lim,"k--","$a = \bar{a}$");
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
rmax = rho-rmin;
nr = 100;
rgrid = linspace(rmin,rmax,nr);


% == Bond Demand == %
S_d = zeros(1,nr);

disp("Calculating net asset supply in the economy, please wait...")
for ir = 1:nr
    [V,con,sav,A] = HJB_solve_new(rgrid(ir));
    [gg,g] = KF_solve_new(A);
    S_d(1,ir) = g(:,1)'*a*da + g(:,2)'*a*da + g(:,3)'*a*da;
end
disp("Done")


% Bond Supply - test1

taxes = tau_0*(y*ydist);

B_d = @(r) taxes./r;

B_dem = B_d(rgrid);


% Bond Supply - test2

% check1 = g(:,1)'*ones(I,1)*da;
% check2 = g(:,2)'*ones(I,1)*da;
% check3 = g(:,3)'*ones(I,1)*da;
% 
% 
% ydist_alt = [check1,check2,check3];
% 
% taxes_alt = tau_0*(y*ydist_alt');
% 
% B_d_alt = @(r) taxes_alt./r;
% 
% B_dem_alt = B_d_alt(rgrid);


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

[B_eq,r_eq] = stationary_equilibrium(ydist);











