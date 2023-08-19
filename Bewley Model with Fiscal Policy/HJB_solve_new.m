%=========================================================================
%                   HJB function solution
%=========================================================================

%=========================================================================
% Solves the system: rho v^{n+1} = u^{n} + A^{n} v^{n+1}
% Returns:
% - value function v
% - Policy functions c,s
% - infinitesimal generator A
%=========================================================================

function [V,con,sav,A] = HJB_solve_new(r)

%=========================================================================
%                   Parameters and calibration
%=========================================================================

% == Household == %

% Poisson Intensities
lambda_12 = 1.5; lambda_13 = 0.5; lambda_21 = 0.75;
lambda_23 = 0.75; lambda_31 = 0.5; lambda_32 = 0.75;

la = [(lambda_12+lambda_13), lambda_12, lambda_13;
      lambda_21, (lambda_21+lambda_23), lambda_23;
      lambda_31, lambda_32, (lambda_31 + lambda_32)]; 

% Income states
y_1 = 0.1; y_2 = 0.15; y_3 = 0.2;

% tax
tau = 0.2;

% borrowing constraint
borrow_lim = -0.15;

% Preferences
gamma = 2; % IES
rho = 0.05; % beta ~ 0.99

%=== technicalities ===%
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

% Iteration criteria
maxit= 100;
crit = 10^(-6);
Delta = 1000;

%=== Building the infinitesimal generator (A Matrix) ===%
% change here
dVf = zeros(I,3); % forward difference approx matrix
dVb = zeros(I,3); % backward difference approx matrix
% c = zeros(I,2); % ?

% Complements the HJB transition matrix with the lambdas
Aswitch = [-speye(I)*la(1,1), speye(I)*la(1,2), speye(I)*la(1,3);
           speye(I)*la(2,1), -speye(I)*la(2,2), speye(I)*la(2,3);
           speye(I)*la(3,1), speye(I)*la(3,2), -speye(I)*la(3,3)];

%INITIAL GUESS - change
v0(:,1) = ((1-tau)*y(1) + r.*a).^(1-gamma)/(1-gamma)/rho; 
v0(:,2) = ((1-tau)*y(2) + r.*a).^(1-gamma)/(1-gamma)/rho;
v0(:,3) = ((1-tau)*y(3) + r.*a).^(1-gamma)/(1-gamma)/rho;

v = v0; % v^{0} that will be an input on the A matrix

for n=1:maxit
    V = v;
    V_n(:,:,n)=V; %update your n-th matrix ~ Changes size every loop

    %====================================================================
    %                   Build Matrix A^{n} 
    %====================================================================

    % == Build the forward difference matrix == %
    dVf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVf(I,:) = ((1-tau).*y + r.*amax).^(-gamma); %will never be used, but impose state constraint a<=amax just in case
    
    % == Build the backward difference matrix == %
    dVb(2:I,:) = (V(2:I,:)-V(1:I-1,:))/da;
    dVb(1,:) = ((1-tau).*y + r.*amin).^(-gamma); %state constraint boundary condition
    
    % Obs: Indicator if function is concave
    % I_concave = dVb > dVf; 
    
    % == Compute policy functions - C and S == % => change 
    %forward difference
    cf = dVf.^(-1/gamma);
    ssf = (1-tau).*yy + r.*aa - cf;
    %backward difference
    cb = dVb.^(-1/gamma);
    ssb = (1-tau).*yy + r.*aa - cb;
    %consumption and derivative of value function at steady state
    c0 = (1-tau).*yy + r.*aa;
    dV0 = c0.^(-gamma);

    % == Upwind Scheme - Indicator functions == %
    If = ssf > 0; %positive drift --> forward difference
    Ib = ssb < 0; %negative drift --> backward difference
    I0 = (1-If-Ib); %at steady state

    % == Upwind Scheme - Finite Difference approx == %
    dV_Upwind = dVf.*If + dVb.*Ib + dV0.*I0; %important to include third term
    c = dV_Upwind.^(-1/gamma);
    u = c.^(1-gamma)/(1-gamma);

    % == Build A^{n} Matrix == %
    X = - min(ssb,0)/da;
    Y = - max(ssf,0)/da + min(ssb,0)/da;
    Z = max(ssf,0)/da;

    % spdiags( elements in the diagonals; diagonals (d= 0 is the main,
    % -1 below; +1 above etc), IxI matrix)
    A1=spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2=spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
    A3=spdiags(Y(:,3),0,I,I)+spdiags(X(2:I,3),-1,I,I)+spdiags([0;Z(1:I-1,3)],1,I,I);

    A = [A1,sparse(I,I),sparse(I,I);
        sparse(I,I),A2,sparse(I,I);
        sparse(I,I),sparse(I,I),A3] + Aswitch;

    % == Obs: check if A satisfies poisson matrix properties == %
    if max(abs(sum(A,2)))>10^(-9)
       disp('Improper Transition Matrix')
       break
    end

    %====================================================================
    %                   Build HJB system of equations
    %====================================================================
    
    % == Build B^{n} Matrix == %
    B = (1/Delta + rho)*speye(3*I) - A;

    % == Build b^{n} vector == %
    u_stacked = [u(:,1);u(:,2);u(:,3)];
    V_stacked = [V(:,1);V(:,2);u(:,3)];
    
    b = u_stacked + V_stacked/Delta;

    % == solve B^{n} v^{n+1} = b^{n} == %
    V_stacked = B\b; 

    % stack solution
    V = [V_stacked(1:I),V_stacked(I+1:2*I),V_stacked((2*I+1):3*I)];

    %====================================================================
    %                   Check for convergence
    %====================================================================
    
    Vchange = V - v;
    v = V;

    dist(n) = max(max(abs(Vchange)));
    if dist(n)<crit
        if display >= 1
            disp('Value Function Converged, Iteration = ')
            disp(n)
        end
        break % FINISH ITERATION IF DIST BELOW THE CONVERGENCE CRITERION
    end
end
con = c;
sav = (1-tau).*yy + r.*aa - con;
end



