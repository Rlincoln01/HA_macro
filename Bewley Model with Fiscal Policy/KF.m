%=========================================================================
%                   Kolmogorov Forward solution
%=========================================================================

%=========================================================================
% Solves the system: rho v^{n+1} = u^{n} + A^{n} v^{n+1}
% Returns:
% - Total Asset Distribution g
% - Marginal densities g_{i} \forall i
%=========================================================================

function [gg,g] = KF(parameters,A)

% display or not message
display = 0;

I = parameters.n_a; %size of the grid
amin = parameters.borrow_lim;
amax = parameters.amax;
da = (amax-amin)/(I-1); % step on asset grid


% === Form the KF system A^{T}g = 0 === %
AT = A'; % differential operator
b = zeros(3*I,1); % vector of zeros

% Obs: dirty fix - impose this to obtain solution in form 
% probability distribution. See moll's HACT appendix page 9
i_fix = 1; %need to fix one value, otherwise matrix is singular
b(i_fix)=.1;
row = [zeros(1,i_fix-1),1,zeros(1,3*I-i_fix)];
AT(i_fix,:) = row;

% === Solve system === %
gg = AT\b;
g_sum = gg'*ones(3*I,1)*da;
gg = gg./g_sum;

% === Obtain marginal distributions === %
g = [gg(1:I),gg(I+1:2*I),gg(2*I+1:3*I)];

check1 = g(:,1)'*ones(I,1)*da;
check2 = g(:,2)'*ones(I,1)*da;
check3 = g(:,3)'*ones(I,1)*da;

checks = [check1,check2,check3];

if display >= 1
    if abs(sum(checks)-1) < 10^(-9)
        disp("Complete and Marginal densities calculated successfully! ")
    else
        disp("Marginal densities don't sum 1 => check for mistakes!")
    end    
end

