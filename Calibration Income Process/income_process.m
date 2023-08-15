%=========================================================================
% Approximation of a Jump-Drift process by a CTMC
%
% Author: Rafael Lincoln 
%=========================================================================

% ==== Parameters ==== %

% rho's 
rho1 = 0;
rho2 = 0;


% Kaplan et al. (2018) calibration:
% Volatility parameters
sigma_t = 1.74; 
sigma_p = 1.53;
% drift parameters
beta_t = 0.761;
beta_p = 0.009;
% poisson arrival rates
lambda_t = 0.08;
lambda_p = 0.007;

% Alternative Calibration: Kaplan et al. (2022) calibration
% 
% z = persistent, u = transitory
% Volatility parameters 
% sigma_p = sqrt(0.239);
% sigma_t = sqrt(1.28);
% 
% % drift parameters
% beta_t = 0.347;
% beta_p = 0.012;
% 
% % poisson arrival rates
% lambda_t = 0.063;
% lambda_p = 0.06;


% Build Grid
tgrid_par = 0.792311684654520;
tgrid_width = 3.78503497352302;

pgrid_width = 5.69800327154487;
pgrid_par = 0.641191531992662;

n_t = 3;
n_p = 11;

t_grid = SymmetricPowerSpacedGrid(n_t,tgrid_par,tgrid_width/2,0);
p_grid = SymmetricPowerSpacedGrid(n_p,pgrid_par,pgrid_width/2,0);
t_steps = diff(t_grid);
p_steps = diff(p_grid);

% ======== 1: Construct Jump Matrices for jump-drift processes ======== %

% 1.1: Construct identity matrices + preallocation

p_eye = eye(n_p); % on Moll's code, l1 = transitory; l2 = persistent
t_eye = eye(n_t);
% preallocation
p_jump = zeros(n_p,n_p);
t_jump = zeros(n_t,n_t);

% Parameters, as in Kaplan et al. (2018)
RestrictJumpDomain = 0;
deltaforapprox = 1;
DriftPointsVersion = 2;

% 1.2: Construct matrix for persistent jump process

for i=1:n_p % row
    for j=1:n_p % column
        if RestrictJumpDomain==0
            if j == 1 % 1st column
                p_jump(i,j) = normcdf(p_grid(j) + (1/2)*p_steps(j),rho1*p_grid(i),sigma_p);
            elseif j>1 && j<n_p % Intermediary columns 
                p_jump(i,j) = normcdf(p_grid(j) + (1/2)*p_steps(j),rho1*p_grid(i),sigma_p) - normcdf(p_grid(j) - (1/2)*p_steps(j-1),rho1*p_grid(i),sigma_p);
            elseif j==n_p % last column
                p_jump(i,j) = 1 - normcdf(p_grid(j)-(1/2)*p_steps(j-1),rho1*p_grid(i),sigma_p);
            end
        elseif RestrictJumpDomain==1
            % define truncated normal
            pd = makedist("Normal","mu",rho1*p_grid(i),"sigma",sigma_p);
            truncated_norm = truncate(pd,p_grid(1),p_grid(n_p));
            if j == 1
                p_jump(i,j) = cdf(truncated_norm,p_grid(j)+(1/2)*p_steps(j));
            elseif j>1 && j < n_p
                p_jump(i,j) = cdf(truncated_norm,p_grid(j)+(1/2)*p_steps(j)) - cdf(truncated_norm,p_grid(j)-(1/2)*p_steps(j-1));
            elseif j==n_p % last column
                p_jump(i,j) = 1 - cdf(truncated_norm,p_grid(j)-(1/2)*p_steps(j-1));
            end    
        end
    end
    p_jump(i,:) = p_jump(i,:)./sum(p_jump(i,:));
end
p_jump = lambda_p.*(p_jump-p_eye);


% 1.3: Construct matrix for transitory jump process

for i=1:n_t % row
    for j=1:n_t % column
        if RestrictJumpDomain==0
            if j == 1 % 1st column
                t_jump(i,j) = normcdf(t_grid(j) + (1/2)*t_steps(j),rho2*t_grid(i),sigma_t);
            elseif j>1 && j<n_t % Intermediary columns 
                t_jump(i,j) = normcdf(t_grid(j) + (1/2)*t_steps(j),rho2*t_grid(i),sigma_t) - normcdf(t_grid(j) - (1/2)*t_steps(j-1),rho1*t_grid(i),sigma_t);
            elseif j==n_t % last column
                t_jump(i,j) = 1 - normcdf(t_grid(j)-(1/2)*t_steps(j-1),rho2*t_grid(i),sigma_t);
            end
        elseif RestrictJumpDomain==1
            % define truncated normal
            pd = makedist("Normal","mu",rho2*t_grid(i),"sigma",sigma_t);
            truncated_norm = truncate(pd,t_grid(1),t_grid(n_t));
            if j == 1
                t_jump(i,j) = cdf(truncated_norm,t_grid(j)+(1/2)*t_steps(j));
            elseif j>1 && j < n_t
                t_jump(i,j) = cdf(truncated_norm,t_grid(j)+(1/2)*t_steps(j)) - cdf(truncated_norm,t_grid(j)-(1/2)*t_steps(j-1));
            elseif j==n_t % last column
                t_jump(i,j) = 1 - cdf(truncated_norm,t_grid(j)-(1/2)*t_steps(j-1));
            end    
        end
    end
    t_jump(i,:) = t_jump(i,:)./sum(t_jump(i,:));
end
t_jump = lambda_t.*(t_jump-t_eye);

% ======== 2: Construct Drift Matrices for the Jump-drift process ======== %

% 2.1 - Build drift matrix for the persistent component
TMdrift_p = zeros(n_p,n_p);

% for i = 1:n_p % rows => Moll's code <-> NOT WORKING
%     if p_grid(i) ~= 0
%         val = (1 - beta_p*deltaforapprox)*p_grid(i);
%         [ii,lp] = FindLinProb1(val,p_grid);
%         if DriftPointsVersion == 2
%             drift_p(i,i) = -1;
%             if p_grid(i) < 0
%                 drift_p(i,i) = drift_p(i,i) + (p_grid(ii(2)) - (1-beta_p)*p_grid(i))/(p_grid(ii(2))-p_grid(i));
%                 drift_p(i,ii(2)) = (-p_grid(i) + (1-beta_p)*p_grid(i))/(p_grid(ii(2))-p_grid(i));
%             elseif p_grid(i) > 0
%                 drift_p(i,ii(1)) = (p_grid(i) - (1-beta_p)*p_grid(i))/(p_grid(i)-p_grid(ii(1)));
%                 drift_p(i,i) = drift_p(i,i) + (-p_grid(ii(1)) + (1-beta_p)*p_grid(i))/(p_grid(i)-p_grid(ii(1)));
%             end    
%         elseif DriftPointsVersion == 3
%             drift_p(i,i) = -1;
%             drift_p(i,ii(1)) = drift_p(i,ii(1)) + (p_grid(ii(2)) - (1-beta_p)*p_grid(i))/(p_grid(ii(2))-p_grid(ii(1)));
%             drift_p(i,ii(2)) = drift_p(i,ii(2)) + (-p_grid(ii(1)) + (1-beta_p)*p_grid(i))/(p_grid(ii(2))-p_grid(ii(1)));
%         end    
%     end
% end


% quick fix for the varying step gri
p_steps(n_p) = p_steps(n_p-1); % this will be useful to calculate the drift matrix

for i = 1:n_p % My code => Provides same result, but works
    drift_p = -beta_p*p_grid(i);
    p_step = p_steps(i);
    up_diag = max(drift_p,0)/p_step;
    low_diag = - min(drift_p,0)/p_step;
    diag = -low_diag - up_diag;
    TMdrift_p(i,i) = diag;
        if i == 1
           TMdrift_p(i,i+1) = up_diag;
        elseif i == n_p
           TMdrift_p(i,i-1) = low_diag;
        else
           TMdrift_p(i,i+1) = up_diag;
           TMdrift_p(i,i-1) = low_diag;
        end
end

% 2.2 - Build drift matrix for the transitory component
TMdrift_t = zeros(n_t,n_t);

% for i = 1:n_t % rows => Moll's code <-> NOT WORKING
%     if t_grid(i) ~= 0
%         val = (1 - beta_t*deltaforapprox)*t_grid(i);
%         [ii,lp] = FindLinProb1(val,t_grid);
%         if DriftPointsVersion == 2
%             drift_t(i,i) = -1;
%             if t_grid(i) < 0
%                 drift_t(i,i) = drift_t(i,i) + (t_grid(ii(2)) - (1-beta_t)*t_grid(i))/(t_grid(ii(2))-t_grid(i));
%                 drift_t(i,ii(2)) = (-t_grid(i) + (1-beta_t)*t_grid(i))/(t_grid(ii(2))-t_grid(i));
%             elseif t_grid(i) > 0
%                 drift_t(i,ii(1)) = (t_grid(i) - (1-beta_t)*t_grid(i))/(t_grid(i)-t_grid(ii(1)));
%                 drift_t(i,i) = drift_t(i,i) + (-t_grid(ii(1)) + (1-beta_t)*t_grid(i))/(t_grid(i)-t_grid(ii(1)));
%             end    
%         elseif DriftPointsVersion == 3
%             drift_t(i,i) = -1;
%             drift_t(i,ii(1)) = drift_t(i,ii(1)) + (t_grid(ii(2)) - (1-beta_t)*t_grid(i))/(t_grid(ii(2))-t_grid(ii(1)));
%             drift_t(i,ii(2)) = drift_t(i,ii(2)) + (-t_grid(ii(1)) + (1-beta_t)*t_grid(i))/(t_grid(ii(2))-t_grid(ii(1)));
%         end    
%     end
% end

for i = 1:n_t % My code => Provides same result, but works
    drift_t = -beta_t*t_grid(i);
    t_step = t_steps(1);
    up_diag = max(drift_t,0)/t_step;
    low_diag = - min(drift_t,0)/t_step;
    diag = -low_diag - up_diag;
    TMdrift_t(i,i) = diag;
        if i == 1
           TMdrift_t(i,i+1) = up_diag;
        elseif i == n_t
           TMdrift_t(i,i-1) = low_diag;
        else
           TMdrift_t(i,i+1) = up_diag;
           TMdrift_t(i,i-1) = low_diag;
        end
end

% Form total Transition Matrix of persistent and transitory shocks

TM_p = p_jump + TMdrift_p;
TM_t = t_jump + TMdrift_t;

% Sanity check: Find ergodic distributions generated by these CTMC

p_dist = stationary_dist(TM_p,n_p);
t_dist = stationary_dist(TM_t,n_t);


% ======== 3: Form a total grid for income process ======== %

% 3.1 - Form the total grid for income process %
% l1 = transitory; l2 = persistent
TM_y = zeros(n_p*n_t,n_p*n_t);
ydist = zeros(n_p*n_t,1);
ygrid = zeros(n_p*n_t,1);

k = 0;
for i=1:n_t
    for j = 1:n_p
        k = k +1;
        ygrid(k) = p_grid(j) + t_grid(i);
        ydist(k) = p_dist(j)*t_dist(i);

        l = 0;
        for m = 1:n_t
            for q = 1:n_p
                l = l+1;
                
%                 ytrans(k,l) = 
                if i == m && j == q
                    TM_y(k,l) = TM_t(i,m) + TM_p(j,q);
                elseif i == m && j ~= q
                    TM_y(k,l) = TM_p(j,q);
                elseif i ~= m && j == q
                    TM_y(k,l) = TM_t(i,m);
                else 
                    TM_y(k,l) = 0;
                end    
            end
        end
    end
end

% 3.2 - Sort in asceding order

ygrid_ordered = sort(ygrid);
for j=1:n_p*n_t
    val = ygrid_ordered(j);
    pos = find(ygrid == val);
    iorder(j) = pos;
end

ygrid_new = zeros(n_p*n_t,1);
TM_y_new = zeros(n_p*n_t,n_p*n_t);
for iy = 1:n_t*n_p
%     ygrid_new(iy) = ygrid(iorder(iy));
    for iy2 = 1:n_p*n_t
        TM_y_new(iy,iy2) = TM_y(iorder(iy),iorder(iy2));
    end
end

TM_y = TM_y_new;
ygrid = ygrid_ordered;

ydist = stationary_dist(TM_y,n_p*n_t);

% ======== 4: Simulate the income process as in Kaplan et al. (2018) ======== %

% Time step (Kaplan's calibration)
dt = 0.25; % time step in quarters
nsim = 5000;
Tsim = floor(20/dt) + 1;

% 4.1 - Discretize Transition Matrix
t_trans = dt*TM_t + t_eye;
p_trans = dt*TM_p + p_eye;

% 4.2 - Define these objects as quarterly MC transition matrices
TM_p_qtr = p_trans;
TM_t_qtr = t_trans;

for it = 1:(floor(1/dt)-1)
    TM_p_qtr = TM_p_qtr*p_trans;
    TM_t_qtr = TM_t_qtr*t_trans; 
end

%cumulative time vector (in quarters)
cumtvec = (1:Tsim) * dt;

% call random numbers for simulation
t_rand = randn([nsim,Tsim]);
p_rand = randn([nsim,Tsim]);

for in=1:nsim
%simulate income path in dt increments
     tsimI(in,1) = gendist(t_dist',1,1);
     psimI(in,1) = gendist(p_dist',1,1);
%     tsimI(in,1) = DiscreteDist1(n_t, t_dist, t_rand(in,1));
%     psimI(in,1) = DiscreteDist1(n_p, p_dist, p_rand(in,1));
     ysim(in,1) = t_grid(tsimI(in,1)) + p_grid(psimI(in,1));

    for it = 2:Tsim
%         tsimI(in,it) = gendist(TM_t_qtr(tsimI(in,it-1),:),1,1);
%         psimI(in,it) = gendist(TM_p_qtr(psimI(in,it-1),:),1,1);
        tsimI(in,it) = gendist(t_trans(tsimI(in,it-1),:),1,1);
        psimI(in,it) = gendist(p_trans(psimI(in,it-1),:),1,1);
%         tsimI(in,it) = DiscreteDist1(n_t, t_trans(tsimI(in,it-1),:), t_rand(in,it));
%         psimI(in,it) = DiscreteDist1(n_p, p_trans(psimI(in,it-1),:), p_rand(in,it));
        ysim(in,it) = t_grid(tsimI(in,it)) + p_grid(psimI(in,it));
    end


    %aggregate to annual income
    ylevsim(in,:) = exp(ysim(in,:));
    for it = 1:5
        yannlevsim(in,it) = sum(ylevsim(in,cumtvec>4.0*(it-1) & cumtvec<=4.0*it));
    end
    yannsim(in,:) = log(yannlevsim(in,:));
end

% ======== 5: Calculate the Moments of the process ======== %

% central moments: logs
muy = sum(yannsim(:,1))/nsim;
mu2y = sum((yannsim(:,1)-muy).^2)/nsim;
mu3y = sum((yannsim(:,1)-muy).^3)/nsim;
mu4y = sum((yannsim(:,1)-muy).^4)/nsim;

% standardised moments: logs
if mu2y > 0
    gam3y = mu3y/(mu2y^1.5);
    gam4y = mu4y/(mu2y^2);
else
    gam3y = 0;
    gam4y = 0;
end

% central moments: 1 year log changes
mudy1 = sum(yannsim(:,2)-yannsim(:,1))/nsim;
mu2dy1 = sum((yannsim(:,2)-yannsim(:,1)-mudy1).^2)/nsim;
mu3dy1 = sum((yannsim(:,2)-yannsim(:,1)-mudy1).^3)/nsim;
mu4dy1 = sum((yannsim(:,2)-yannsim(:,1)-mudy1).^4)/nsim;

kurtdy1 = mu4dy1/(mu2dy1^2);

% standardised moments: 1 year log changes
if mu2dy1 > 0
gam3dy1 = mu3dy1/(mu2dy1^1.5);
gam4dy1 = mu4dy1/(mu2dy1^2);
else
gam3dy1 = 0;
gam4dy1 = 0;
end

% central moments: 5 year log changes
mudy5 = sum(yannsim(:,5)-yannsim(:,1))/nsim;
mu2dy5 = sum((yannsim(:,5)-yannsim(:,1)-mudy5).^2)/nsim;
mu3dy5 = sum((yannsim(:,5)-yannsim(:,1)-mudy5).^3)/nsim;
mu4dy5 = sum((yannsim(:,5)-yannsim(:,1)-mudy5).^4)/nsim;

kurtdy5 = mu4dy5/(mu2dy5^2);

% standardised moments: 5 year log changes
if mu2dy5 > 0
gam3dy5 = mu3dy5/(mu2dy5^1.5);
gam4dy5 = mu4dy5/(mu2dy5^2);
else
gam3dy5 = 0;
gam4dy5 = 0;
end

% fraction 1 year log changes in ranges
fracdy1less5 = sum(abs(yannsim(:,2)-yannsim(:,1)) < 0.05)/double(nsim);
fracdy1less10 = sum(abs(yannsim(:,2)-yannsim(:,1)) < 0.1)/double(nsim);
fracdy1less20 = sum(abs(yannsim(:,2)-yannsim(:,1)) < 0.2)/double(nsim);
fracdy1less50 = sum(abs(yannsim(:,2)-yannsim(:,1)) < 0.5)/double(nsim);

% ======== 6: Make plots ======== %

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
sigmady1 = sqrt(mu2dy1);
sigmady5 = sqrt(mu2dy5);
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
st_dev_text = sprintf('St. Dev = %0.2f', mu2dy1);
kurtosis_text = sprintf('Kurtosis = %0.2f', kurtdy1);
dim = [0.15 0.15 0.2 0.2];
annotation('textbox',dim,'String',{st_dev_text, kurtosis_text},'FitBoxToText','on', ...
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
st_dev_text = sprintf('St. Dev = %0.2f', mu2dy5);
kurtosis_text = sprintf('Kurtosis = %0.2f', kurtdy5);
dim = [0.6 0.15 0.2 0.2];
annotation('textbox',dim,'String',{st_dev_text, kurtosis_text},'FitBoxToText','on', ...
    "Interpreter","latex","FontSize",14,'BackgroundColor','white');




% ======== 5: Simulate the income process -> My code ======== %


% [jmptout,sttout] = simCTMC(TM_y,1000,5000);
% 
% times = linspace(0,1000,1000);
% 
% idxs = simCTMC_t(times,jmptout,sttout);
% 
% y_series = ygrid(idxs);
% 
% dy_series = diff(y_series,1);
% 
% % mean = moment(y_series, 1, "all"); %mean
% % variance = moment(y_series, 2, "all"); %variance
% % skew = moment(y_series, 3, "all"); % skewness
% % kurtosis = moment(y_series, 4, "all"); % kurtosis
% 
% % dy_5_series = diff(y_series,1,2);
% 
% [f, xi] = ksdensity(dy_series(:),ygrid);
% % plot a normal distribution with variance of the data to compare
% y_norm = normpdf(xi,0,variance);
% 
% 
% fig= figure;
% fig.Position = [151,208,700,585];
% set(gca,'FontSize',16);
% p1 = plot(xi, f, "LineWidth",2);
% p1.Color = '#000080';
% hold on
% p2 = plot(xi,y_norm,"LineStyle","--","LineWidth",2);
% p2.Color = "r";
% lgd = legend(["Generated by Income distribution","Generated by $\mathcal{N}(0,\sigma^{2})$"],...
%     AutoUpdate="off",Location="northeast");
% lgd.Interpreter = "latex";
% lgd.FontSize = 14;
% tlt = title("$\ln y_{t+1} - \ln y_{t}$","FontSize",20);
% tlt.Interpreter = "latex";
% 
% % Test: Import their markov matrix
% 
% filename = 'ymarkov_combined.txt';
% TM = dlmread(filename);
% 
% [jmptout,sttout] = simCTMC(TM,1000,5000);
% 
% times = linspace(0,1000,1000);
% 
% idxs = simCTMC_t(times,jmptout,sttout);
% 
% y_series = ygrid(idxs);
% 
% dy_series = diff(y_series,1);
% 
% mean = moment(y_series, 1, "all"); %mean
% variance = moment(dy_series, 2, "all"); %variance
% skew = moment(y_series, 3, "all"); % skewness
% kurtosis = moment(y_series, 4, "all"); % kurtosis
% 
% % dy_5_series = diff(y_series,1,2);
% 
% [f, xi] = ksdensity(dy_series(:),ygrid);
% % plot a normal distribution with variance of the data to compare
% y_norm = normpdf(xi,0,variance);
% 
% 
% fig= figure;
% fig.Position = [151,208,700,585];
% set(gca,'FontSize',16);
% p1 = plot(xi, f, "LineWidth",2);
% p1.Color = '#000080';
% hold on
% p2 = plot(xi,y_norm,"LineStyle","--","LineWidth",2);
% p2.Color = "r";
% lgd = legend(["Generated by Income distribution","Generated by $\mathcal{N}(0,\sigma^{2})$"],...
%     AutoUpdate="off",Location="northeast");
% lgd.Interpreter = "latex";
% lgd.FontSize = 14;
% tlt = title("$\ln y_{t+1} - \ln y_{t}$","FontSize",20);
% tlt.Interpreter = "latex";
% 
% 
% 
% 
% 


