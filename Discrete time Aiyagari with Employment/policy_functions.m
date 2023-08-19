%=========================================================================
%                     Consumption and saving policy functions
%=========================================================================

function [con,sav] = policy_functions(r,w,na,borrow_lim,amax) 
    %% Set Parameters
    
    % Household Parameters
    beta = 0.99;
    
    % Capital Accumulation
    delta = 0.025;
    
    % return on asset
    R = 1 + r -delta; %gross return on capital
    
    % Labor Market
    l_bar = 1/0.9; % # of hours worked
    tau = 0.015; % tax on earnings
    mu_unemp = 0.15; % share of wage received as unemployment benefit
    
    % preferences
    u = @(c)log(c); %utility function
    
    u1 = @(c)1./c; % first derivative of the utility function
    
    u1inv = @(u) 1./u; % inverse of the first derivative of the utility function

    
    
    %% Set up  grids
    
    % asset grid
    agrid_par   = 0.4; %1 for linear, 0 for L-shaped
    agrid = linspace(0,1,na)';
    agrid = agrid.^(1./agrid_par);
    agrid = borrow_lim + (amax-borrow_lim).*agrid;
    
    % income: markov chain with employment
    p_00 = 0.6; % probability of remaining unemployed
    p_01 = 1-p_00; % probability of being employed
    p_11 = 0.955555; % probability of remaining employed
    p_10 = 1 - p_11; % probability of being fired
    
    P = [p_00 p_01; p_10 p_11]; % transition matrix
    
    % Employment grid
    ne = 2; % n of points in the grid
    y = @(e) e*(1-tau)*l_bar*w + (1-e)*mu_unemp*w; % returns labor income given if unemployed or no
    
   
    
    %% Computation parametrization
    
    % computation
    max_iter    = 1000;
    tol_iter    = 1.0e-6;
    iter = 0;
    cdiff = 1000;
    Display     = 0;
    
    
    %% Initialize consumption function
    
    conguess = zeros(na,ne);
    for ie = 1:ne
        conguess(:,ie) = (r-delta).*agrid + y(ie-1);
    end
    
    %% Iterate on Euler Equation with the endogenous grid point methods
    
    con = conguess;
    
    while iter <= max_iter && cdiff>tol_iter
        iter = iter + 1;
        conlast = con;
        sav = zeros(na,ne);
        ass1 = zeros(na,ne);
        for ie = 1:ne %fix income first in order to establish endogenous asset grid
            % determine the income distribution through a markov chain
            if ie == 1 % if he's unemployed
                    p = [1 0];
            elseif ie == 2 % if he's employed
                    p = [0 1];
            end
                ydist = transpose(P)*transpose(p); %conditional distribution
                emuc = u1(con)*ydist;
                muc1 = beta.*R.*emuc; %expected marginal utility of consumption
                con1 = u1inv(muc1); % solving for the consumption policy function (from y and a tomorrow)
                ass1(:,ie) = (con1 + agrid -y(ie-1))./R; % asset today as a function of income and assets tomorrow
                for ia = 1:na % Partition the grid into the point above and below the optimal asset today
                    if agrid(ia) < ass1(1,ie) % asset today evaluated at the future borrowing constraint
                        sav(ia,ie) = borrow_lim; % savings = assets = borroming limit
                    else % if not, assets today are a linear interpolation of assets tomorrow and gridpoint
                        sav(ia,ie) = lininterp1(ass1(:,ie),agrid,agrid(ia));
                    end
                end
            con(:,ie) = R.*agrid +y(ie-1) - sav(:,ie);
        end
    
        cdiff = max(max(abs(con-conlast)));
        if Display >= 1 && mod(iter,20) ==0
            disp(['Iteration no. ' int2str(iter), ' max con fn diff is ' num2str(cdiff)]);
        end
    end
end





