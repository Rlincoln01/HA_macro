% ========================================================================
%                   Loss Function
% ========================================================================
% Description: 
%
% 
%
%
% Author: Rafael Lincoln
% ========================================================================


function [loss,decomp] = loss_fct(sample_moments,sim_moments,weights)

% Take the percentage deviation of the empirical and simulated moment
moment_names = fieldnames(sim_moments);
n_f = length(moment_names);
F = zeros(n_f,1);


% check if weights are provided separately or not
if nargin < 3
    % if not provided, I use an identity matrix
    weights = ones(1,n_f);
end


for i =1:n_f
    moment = moment_names{i};
    % percent deviation from the sample moment
    F(i) = (sim_moments.(moment) - sample_moments.(moment))/abs(sample_moments.(moment));
    % Absolute deviation from the sample moment
    % F(i) = (sim_moments.(moment) - sample_moments.(moment));
end

% weighting matrix
W = weights.*eye(n_f);

% loss function
loss = F'*W*F;

% decomposition of deviations
dec = 100.*(weights'.*(F.^2))./loss;

for i =1:n_f
    moment = moment_names{i};
    % percent deviation from the sample moment
    decomp.(moment) = dec(i);
end

end










