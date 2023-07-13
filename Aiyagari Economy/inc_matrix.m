% ========================================================================
%                   Income Process Generator Matrix
% ========================================================================
% Description: Generates the CTMC transition matrix for the jump-drift 
% process of productivity of households
%
% permanent_params = [n_p,beta_p,sigma_p,lambda_p]
% trans_params = [n_t,beta_t,sigma_t,lambda_t]
%
% Author: Rafael Lincoln / Kaplan et al. (2018) (adapted from their fortran
% code)
%
% ========================================================================

function [TM_y,ygrid,TM_p,p_grid,TM_t,t_grid] = inc_matrix(beta_p,beta_t,sigma_p,sigma_t,lambda_p,lambda_t)

global n_p n_t pgrid_width pgrid_par tgrid_width tgrid_par

% ========== Parameters ========== %

% rho's => is the shock centered at that gridpoint or no?
rho1 = 0;
rho2 = 0;


% ========== Grid ========== %

t_grid = SymmetricPowerSpacedGrid(n_t,tgrid_par,tgrid_width,0);
p_grid = SymmetricPowerSpacedGrid(n_p,pgrid_par,pgrid_width,0);

% Measure of grid
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

% quick fix for the varying step gri
p_steps(n_p) = p_steps(n_p-1); % this will be useful to calculate the drift matrix

for i = 1:n_p 
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

for i = 1:n_t 
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
% 
% p_dist = stationary_dist(TM_p,n_p);
% t_dist = stationary_dist(TM_t,n_t);

% ======== 3: Form a total grid for income process ======== %

% 3.1 - Form the total grid for income process %
% l1 = transitory; l2 = persistent
TM_y = zeros(n_p*n_t,n_p*n_t);
% ydist = zeros(n_p*n_t,1);
ygrid = zeros(n_p*n_t,1);

k = 0;
for i=1:n_t
    for j = 1:n_p
        k = k +1;
        ygrid(k) = p_grid(j) + t_grid(i);
        % ydist(k) = p_dist(j)*t_dist(i);
        l = 0;
        for m = 1:n_t
            for q = 1:n_p
                l = l+1;
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


TM_y_new = zeros(n_p*n_t,n_p*n_t);
for iy = 1:n_t*n_p
%     ygrid_new(iy) = ygrid(iorder(iy));
    for iy2 = 1:n_p*n_t
        TM_y_new(iy,iy2) = TM_y(iorder(iy),iorder(iy2));
    end
end

TM_y = TM_y_new;
ygrid = ygrid_ordered;
% ydist = stationary_dist(TM_y,n_p*n_t);

end

%=========================================================================
% SymmetricPowerSpacedGrid
%
% Gives a grid spaced between center-width and center+width based on the 
% interval [-1,1] with a function x^(1/k) either side. If k = 1, the grid
% is linear; and if k = 0, is L-shaped
%
% Obs: n>= 2
% Author: Greg Kaplan, G. Violante and Ben Moll
%=========================================================================

function[y] = SymmetricPowerSpacedGrid(n,k,width,center)

if n == 2
    y(1) = center - width;
    y(2) = center + width;
end

% create a linear spaced grid between -1,1
x = linspace(-1,1,n);

for i=1:n
    if x(i)>0
        z(i) = x(i)^(1/k);
    elseif x(i) == 0
        z(i) = 0;
    else 
        z(i) = -((-x(i))^(1/k));
    end
y= center + width.*z;    
end
end