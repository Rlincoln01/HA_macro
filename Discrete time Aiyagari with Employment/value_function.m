%=========================================================================
%                     Value function
%=========================================================================

function [V] = value_function(con,cap,na,ne,P,agrid)
% == Value Function == %
u = @(c)log(c); %utility function

% params
beta = 0.99;

Vguess = u(con)./(1-beta);

max_iter    = 1000;
tol_iter    = 1.0e-6;
    
iter = 0;
Vdiff = 1000;


V = Vguess;
while iter < max_iter && Vdiff>tol_iter    
    V_next = zeros(na,ne);
    for ie = 1:ne
        % Iterate value function
        V_next(:,ie) = u(con(:,ie))+ beta*P(ie,1)*interp1(agrid,V(:,1),cap(:,ie),"linear","extrap")...
            + beta*P(ie,2)*interp1(agrid,V(:,2),cap(:,ie),"linear","extrap");
    end
    Vdiff = max(max(abs(V_next-V)));
    V = V_next;
    iter = iter +1;    
end
