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


function loss = loss_fct(sample_moments,sim_moments)
% Take the percentage deviation of the empirical and simulated moment
moment_names = fieldnames(sim_moments);
n_f = length(moment_names);
F = zeros(n_f,1);

for i =1:n_f
    moment = moment_names{i};
    % percent deviation from the sample moment
    F(i) = (sim_moments.(moment) - sample_moments.(moment))/abs(sample_moments.(moment));
end

% weighting matrix
W = eye(n_f);

% loss function
loss = F'*W*F;

if isnan(loss)
    loss = 10^6;
end

end










