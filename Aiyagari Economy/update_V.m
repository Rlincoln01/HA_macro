%=========================================================================
%                   Update value function with the implicit method
%=========================================================================
% Description: Discretization of the Infinitesimal generator with the
%              upwind method
% 
%=========================================================================


function V = update_V(V_initial,u_n,L,n_z,n_a,delta)

global rho zeta

% =========================== Update V ========================== %

% 4.1 - Build B matrix
B = (1+delta*(rho+ zeta))*speye(n_a*n_z) - delta*L;

% 4.2 - Build b^n vector
u_n_stacked = reshape(u_n,n_a*n_z,1);
V_stacked = reshape(V_initial,n_a*n_z,1);
b = delta*u_n_stacked + V_stacked;

% solve for V^{n+1}
V_stacked = B\b;

% reestack solution
V = reshape(V_stacked,n_a,n_z);
end