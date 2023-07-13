%=========================================================================
%                   Discretize the Infinitesimal generator
%=========================================================================
% Description: Discretization of the Infinitesimal generator with the
%              upwind method
% 
%=========================================================================

function [A,c] = Build_A(V,w,r,lump_sum,omega,a_grid,z_grid)

global amax tau gamma

% ======================= Build the grid ======================= %

% Number of points
n_a = length(a_grid);
n_z = length(z_grid);
amin = a_grid(1); % borrowing constraint

% Build varying measure for the asset grid
[~, daf, dab] = grid_measure(a_grid,n_a);

daaf = daf*ones(1,n_z);           % Matrix IXJ of F-measures
daab = dab*ones(1,n_z);           % Matrix IXJ of B-measures

% 4.1) Determine joint grid of assets and income
zz = ones(n_a,1)*z_grid';                   % Matrix IxM of productivity
aa = a_grid*ones(1,n_z);                  % Matrix IxM of assets

% ======================= Build the A matrix ==================== %

% 1.1 - Build the Forward difference matrix

Vaf(1:n_a-1,:) = max(exp(1)^(-10),(V(2:n_a,:)-V(1:n_a-1,:))./daaf(1:n_a-1,:));
Vaf(n_a,:) = ((1-tau)*w*exp(z_grid + omega) + r.*amax + lump_sum).^(-gamma);

% 1.2 - Build the Backward difference matrix

Vab(2:n_a,:) = max(exp(1)^(-10),(V(2:n_a,:)-V(1:n_a-1,:))./daab(2:n_a,:));
Vab(1,:) = (max(exp(1)^(-10),(1-tau)*w*exp(z_grid + omega) + r.*amin + lump_sum)).^(-gamma);  % state constraint boundary condition

% 1.3 - Build F&B Policy Functions

% forward difference policy functions
cf = max(exp(1)^(-10),Vaf.^(-1/gamma));
sf = (1-tau)*w*exp(zz + omega) + r.*aa + lump_sum- cf;

% backward difference policy functions
cb = max(exp(1)^(-10),Vab.^(-1/gamma));
sb = (1-tau)*w*exp(zz + omega) + r.*aa + lump_sum - cb;

% consumption and derivative of value function at steady state
c0 = max(exp(1)^(-10),(1-tau)*w*exp(zz + omega) + r.*aa + lump_sum);
Va0 = c0.^(-gamma);

% 1.4 - Upwind Scheme

% Indicator Functions
If = sf > 0; % positive drift --> forward difference
Ib = sb < 0; % negative drift --> backward difference
I0 = (1-If-Ib); % at steady state

% 1.5 - Building the Finite Difference approximation
Va_Upwind = Vaf.*If + Vab.*Ib + Va0.*I0; %important to include third term

% 1.6 - Building FD aooriximation of Policy
c = max(exp(1)^(-10),Va_Upwind.^(-1/gamma));


% 1.7 - Build the sparse A^n matrix
X = - min(sb,0)./daab;
Y = - max(sf,0)./daaf + min(sb,0)./daab;
Z = max(sf,0)./daaf;

% 1.7.1 - Upper diagonal
updiag=0; %This is needed because of the peculiarity of spdiags.
for m=1:n_z
    updiag=[updiag;Z(1:n_a-1,m);0];
end

% 1.7.2 - Build central diagonal
centdiag=reshape(Y,n_a*n_z,1);

% 1.7.3 - Build lower diagonal
lowdiag=X(2:n_a,1);
for m=2:n_z
    lowdiag=[lowdiag;0;X(2:n_a,m)];
end

% 1.7.4 - Assemble diagonals as a sparse matrix
A =spdiags(centdiag,0,n_a*n_z,n_a*n_z)+spdiags([updiag;0],1,n_a*n_z,n_a*n_z)+spdiags([lowdiag;0],-1,n_a*n_z,n_a*n_z);

% == 3rd) Check if it is a transition matrix (sum(row) = 1, for all rows)
if max(abs(sum(A,2)))>10^(-9)
    disp('Improper Transition Matrix')
end

end    