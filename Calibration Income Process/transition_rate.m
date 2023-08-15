function r = transition_rate(x, y, dt, beta_xi, sigma_eta)
    a = (x - beta_xi*y)*dt;
    b = sqrt(dt)*sigma_eta;
    r = (1/sqrt(2*pi*b^2))*exp(-(y - exp(-a)*x)^2/(2*b^2));
end