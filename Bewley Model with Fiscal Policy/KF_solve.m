%=========================================================================
%                   Kolmogorov Forward solution
%=========================================================================

%=========================================================================
% Solves the system: rho v^{n+1} = u^{n} + A^{n} v^{n+1}
% Returns:
% - Total Asset Distribution g
% - Marginal densities g_{i} \forall i
%=========================================================================

function [gg,g] = KF_solve(A)

I=500; % size of asset grid
borrow_lim = -0.15;
amin = borrow_lim;
amax = 3;
da = (amax-amin)/(I-1); % step on asset grid


% === Form the KF system A^{T}g = 0 === %
AT = A'; % differential operator
b = zeros(2*I,1); % vector of zeros

% Obs: dirty fix - impose this to obtain solution in form 
% probability distribution. See moll's HACT page 9
i_fix = 1; %need to fix one value, otherwise matrix is singular
b(i_fix)=.1;
row = [zeros(1,i_fix-1),1,zeros(1,2*I-i_fix)];
AT(i_fix,:) = row;

% === Solve system === %
gg = AT\b;
g_sum = gg'*ones(2*I,1)*da;
gg = gg./g_sum;

% === Obtain marginal distributions === %
g = [gg(1:I),gg(I+1:2*I)];
