% ========================================================================
%                   Moments
% ========================================================================
% Description: calculates the moments of a set of simulated data
%              
%
% Inputs:
%   - Yearly income in log for each income [Matrix (n_sim,t_sim)]
%
%
% Author: Rafael Lincoln 
% ========================================================================

function sim_moments = moments(yannsim,settings)

nsim = settings.n_sim;          % Number of income path simulations
t=0;

% central moments: logs
muy = sum(yannsim(:,1+t))/nsim;
mu2y = sum((yannsim(:,1+t)-muy).^2)/nsim;
mu3y = sum((yannsim(:,1+t)-muy).^3)/nsim;
mu4y = sum((yannsim(:,1+t)-muy).^4)/nsim;

% standardised moments: logs
if mu2y > 0
    gam3y = mu3y/(mu2y^1.5);
    gam4y = mu4y/(mu2y^2);
else
    gam3y = 0;
    gam4y = 0;
end

% central moments: 1 year log changes
mudy1 = sum(yannsim(:,t+2)-yannsim(:,t+1))/nsim;
mu2dy1 = sum((yannsim(:,t+2)-yannsim(:,t+1)-mudy1).^2)/nsim;
mu3dy1 = sum((yannsim(:,t+2)-yannsim(:,t+1)-mudy1).^3)/nsim;
mu4dy1 = sum((yannsim(:,t+2)-yannsim(:,t+1)-mudy1).^4)/nsim;

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
mudy5 = sum(yannsim(:,t+5)-yannsim(:,t+1))/nsim;
mu2dy5 = sum((yannsim(:,t+5)-yannsim(:,t+1)-mudy5).^2)/nsim;
mu3dy5 = sum((yannsim(:,t+5)-yannsim(:,t+1)-mudy5).^3)/nsim;
mu4dy5 = sum((yannsim(:,t+5)-yannsim(:,t+1)-mudy5).^4)/nsim;

kurtdy5 = mu4dy5/(mu2dy5^2);

% standardised moments: 5 year log changes
if mu2dy5 > 0
gam3dy5 = mu3dy5/(mu2dy5^1.5);
gam4dy5 = mu4dy5/(mu2dy5^2);
else
gam3dy5 = 0;
gam4dy5 = 0;
end

% P5010 and P9050 of the 1 and 5 year log change
pct_1yr_change = prctile(yannsim(:,t+2)-yannsim(:,t+1),[10,50,90]);
pct_5yr_change = prctile(yannsim(:,t+5)-yannsim(:,t+1),[10,50,90]);

% P9050_1yr_change = pct_1yr_change(3) - pct_1yr_change(2);
% P9050_5yr_change = pct_5yr_change(3) - pct_5yr_change(2);
% P5010_1yr_change = pct_1yr_change(2) - pct_1yr_change(1);
% P5010_5yr_change = pct_5yr_change(2) - pct_5yr_change(1);

P9010_1yr_change = pct_1yr_change(3) - pct_1yr_change(1);
P9010_5yr_change = pct_5yr_change(3) - pct_5yr_change(1);


% fraction 1 year log changes in ranges
% fracdy1less5 = sum(abs(yannsim(:,2)-yannsim(:,1)) < 0.05)/double(nsim);
% fracdy1less10 = sum(abs(yannsim(:,2)-yannsim(:,1)) < 0.1)/double(nsim);
% fracdy1less20 = sum(abs(yannsim(:,2)-yannsim(:,1)) < 0.2)/double(nsim);
% fracdy1less50 = sum(abs(yannsim(:,2)-yannsim(:,1)) < 0.5)/double(nsim);

% return desired moments

sim_moments.std_log_y = sqrt(mu2y);
sim_moments.std_1yr_log_change = sqrt(mu2dy1);
sim_moments.std_5yr_log_change = sqrt(mu2dy5);
sim_moments.kurt_1yr_log_change = kurtdy1;
sim_moments.kurt_5yr_log_change = kurtdy5;
% sim_moments.P9050_1yr_change = P9050_1yr_change;
% sim_moments.P9050_5yr_change = P9050_5yr_change;
% sim_moments.P5010_1yr_change = P5010_1yr_change;
% sim_moments.P5010_5yr_change = P5010_5yr_change;
sim_moments.P9010_1yr_log_change = P9010_1yr_change;
sim_moments.P9010_5yr_log_change = P9010_5yr_change;
end
